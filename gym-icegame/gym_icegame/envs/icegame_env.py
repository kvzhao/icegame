from __future__ import division
import gym
from gym import error, spaces, utils, core

#from six import StingIO
import sys, os
import json
import random
import numpy as np
from icegame import SQIceGame, INFO

import time

rnum = np.random.randint

class IceGameEnv(core.Env):
    def __init__ (self, L, kT, J):
        self.L = L
        self.kT = kT
        self.J = J
        self.N = L**2
        num_neighbors = 2
        num_replicas = 1
        num_mcsteps = 2000
        num_bins = 1
        num_thermalization = num_mcsteps
        tempering_period = 1

        self.mc_info = INFO(self.L, self.N, num_neighbors, num_replicas, \
                num_bins, num_mcsteps, tempering_period, num_thermalization)

        self.sim = SQIceGame(self.mc_info)
        self.sim.set_temperature (self.kT)
        self.sim.init_model()
        self.sim.mc_run(num_mcsteps)

        self.episode_terminate = False
        self.accepted_episode = False

        self.name_mapping = dict({
                                  0 :   'right',
                                  1 :   'down',
                                  2 :   'left',
                                  3 :   'up',
                                  4 :   'lower_next',
                                  5 :   'upper_next',
                                  6 :   'metropolis',
                                  })

        self.index_mapping = dict({
                                  'right': 0,
                                  'down' : 1,
                                  'left' : 2,
                                  'up' : 3,
                                  'lower_next' : 4,
                                  'upper_next' : 5,
                                  'metropolis' : 6,
                                  })


        ### action space and state space
        self.action_space = spaces.Discrete(len(self.name_mapping))
        self.observation_space = spaces.Box(low=-1.0, high=1.0, shape=(self.L, self.L, 4))
        self.reward_range = (-1, 1)

        # output file
        self.ofilename = 'loop_sites.log'
        # render file
        self.rfilename = 'loop_renders.log'
        # save log to json for future analysis
        self.json_file = 'env_history.json'

        self.stacked_axis = 2

        ## counts reset()
        self.episode_counter = 0

        ## ray test, add list:
        #     1. log 2D (x, y) in self.ofilename
        #     2. add self.calculate_area() and loop_area
        #     3. auto_metropolis (uncompleted)
        #   8/2:
        #     5. add save_record_dict(): to record the amount of length and area
        #     6. add area_reward in step()
        #     7. add hundred_test(): to count accepted ratio in last 100 episodes
        #     8. add now_position & path flag in get_obs()
        self.record_dict = [{}, {}]
        self.accepted_in_hundred = 0
        self.accepted_in_hundred_stack = []

        ## ray test, step_setting
        self.area_reward = True                 # use both length & area to calculate the reward
        self.auto_metropolis = False            # if the condition is OK, auto execute metropolis
        self.metropolis_terminal = False        # if metropolis_executed == True, terminal = True
        self.strict_step = False                # if rets[0] (is_aceept) == -1, terminal = True


    def step(self, action):
        terminate = False
        reward = 0.0
        obs = None
        rets = [0.0, 0.0, 0.0, 0.0]
        metropolis_executed = False
        self.hundred_test()     ## ray test

        ## execute different type of actions
        if (action == 6):
            self.sim.flip_trajectory()
            rets = self.sim.metropolis()
            metropolis_executed = True
            if (self.metropolis_terminal):
                terminate = True
        elif (0 <= action < 6) :
            rets = self.sim.draw(action)
        
        is_aceept, dEnergy, dDensity, dConfig = rets

        if (self.auto_metropolis):
            if is_aceept > 0 and dConfig > 0:
                metropolis_executed = True

        # metropolis judgement
        if (metropolis_executed):
            if is_aceept > 0 and dConfig > 0:
                self.sim.update_config()
                print ('[GAME_ENV] PROPOSAL ACCEPTED!')
                total_steps = self.sim.get_total_steps()
                ep_steps = self.sim.get_ep_step_counter()
                ep = self.sim.get_episode()
                loop_length = self.sim.get_accepted_length()[-1]
                loop_area = self.calculate_area()
                self.save_record_dict(loop_length, loop_area)   ## ray test
                update_times = self.sim.get_updated_counter()

                reward = 1.0 * (loop_length / 4.0) # reward with different length by normalizing with len 4 elements

                ## ray test, area reward & hundred_test()
                self.hundred_test(True)     ## acceptance ratio in hundred Eps
                if self.area_reward:        ## add area_reward
                    reward = reward + loop_area

                # output to self.ofilename
                with open(self.ofilename, 'a') as f:
                    f.write('1D: {}  (2D: {})\n'.format(self.sim.get_trajectory(), self.conver_1Dto2D(self.sim.get_trajectory())))
                    print ('\tSave loop configuration to file: {}'.format(self.ofilename))

                print ('\tTotal accepted number = {}'.format(self.sim.get_updated_counter()))
                print ('\tAccepted loop length = {}, area = {}'.format(loop_length, loop_area))
                print ('\tAccepted record: loop = {}, area = {}'.format(self.record_dict[0], self.record_dict[1]))  ## ray test
                print ('\tAgent walks {} steps in episode, action counters: {}'.format(ep_steps, self.sim.get_ep_action_counters()))
                action_counters = self.sim.get_action_statistics()
                action_stats = [x / total_steps for x in action_counters]
                print ('\tStatistics of actions all episodes (ep={}, steps={}): \n\t {}'.format(ep, total_steps, [format(member, '.3f') for member in action_stats]))
                print ('\tAcceptance ratio (accepted/total Eps) = {} %, Acceptance ratio in hundred = {} %'.format(update_times * 100.0 / ep, self.accepted_in_hundred))
                print ('-' + self.L * '---' + '+\n')
                self.dump_env_states()

                self.render(info=[loop_length, loop_area])
                self.sim.clear_buffer()
            else:
                self.sim.clear_buffer()
                terminate = True
                # Avoid running metropolis at start
                if (rets[3] == 0.0):
                    reward = -0.8
            # reset or update
        else:
            reward = self._stepwise_weighted_returns(rets)
            # as usual
            if (self.strict_step): 
                if is_aceept < 0:
                    terminate = True

        obs = self.get_obs()
        ## add timeout mechanism?

        return obs, reward, terminate, rets

    def _stepwise_weighted_returns(self, rets):
        icemove_w = 0.000
        energy_w = -1.0
        defect_w = 0.0
        baseline = 0.009765625 ## 1 / 1024
        scaling = 2.0
        return (icemove_w * rets[0] + energy_w * rets[1] + defect_w * rets[2] + baseline) * scaling

    # Start function used for agent learing
    def start(self, init_site=None):
        if init_site == None:
            init_agent_site = self.sim.start(rnum(self.N))
        else:
            init_agent_site = self.sim.start(init_site)
        assert(init_site == init_agent_site)

    def reset(self):
        ## clear buffer and set new start of agent
        site = rnum(self.N)
        init_site = self.sim.restart(site)
        assert(init_site == site)
        self.episode_counter += 1
        return self.get_obs()

    def timeout(self):
        return self.sim.timeout()

    @property
    def agent_site(self):
        return self.sim.get_agent_site()

    @property
    def action_name_mapping(self):
        return self.name_mapping

    @property
    def name_action_mapping(self):
        return self.index_mapping


    def render(self, info, mapname ='traj', mode='ansi', close=False):
        #of = StringIO() if mode == 'ansi' else sys.stdout
        #print ('Energy: {}, Defect: {}'.format(self.sqice.cal_energy_diff(), self.sqice.cal_defect_density()))
        s = None
        if (mapname == 'traj'):
            s = self._transf2d(self.sim.get_canvas_map())
        start = self.sim.get_start_point()
        start = (int(start/self.L), int(start%self.L))
        s[start] = 3
        screen = '\r'
        screen += '\n\t'
        screen += '+' + self.L * '---' + '+\n'
        for i in range(self.L):
            screen += '\t|'
            for j in range(self.L):
                p = (i, j)
                spin = s[p]
                if spin == -1:
                    screen += ' o '
                elif spin == +1:
                    screen += ' * '
                elif spin == 0:
                    screen += '   '
                elif spin == +2:
                    screen += ' @ '
                elif spin == -2:
                    screen += ' O '
                elif spin == +3:
                    screen += ' x '
            screen += '|\n'
        screen += '\t+' + self.L * '---' + '+\n'
        #sys.stdout.write(screen)
        with open(self.rfilename, 'a') as f:
            f.write('Episode: {}, global step = {}\n'.format(self.episode_counter, self.sim.get_total_steps()))
            f.write('Accepted loop length = {}, area = {}'.format(info[0], info[1]))    ## add
            f.write('{}\n'.format(screen))


    def get_obs(self):  ## use 1D coordinate

        ## ray test, now_position flag in 'canvas_map'
        # in '_canvas_map':
        #   now_position -> now_position_flag
        #   path -> path_flag

        now_position_flag = +1.0    ## try -1
        path_flag = -1.0
        now_position = self.sim.get_agent_site()
        _canvas_map = np.copy(self.sim.get_canvas_map())
        for index in np.where(_canvas_map==now_position_flag)[0]:
            _canvas_map[index] = path_flag

        if _canvas_map[self.sim.get_start_point()] == 0:
            _canvas_map[self.sim.get_start_point()] = path_flag

        _canvas_map[now_position] = now_position_flag

        config_map = self._transf2d(self.sim.get_state_t_map())
        canvas_map = self._transf2d(_canvas_map)
        energy_map = self._transf2d(self.sim.get_energy_map())
        defect_map = self._transf2d(self.sim.get_defect_map())

        return np.stack([config_map,
                         canvas_map,
                         energy_map,
                         defect_map
        ], axis=self.stacked_axis)

    @property
    def unwrapped(self):
        """Completely unwrap this env.
            Returns:
                gym.Env: The base non-wrapped gym.Env instance
        """
        return self

    def sysinfo(self):
        print ('')

    def _transf2d(self, s):
        # do we need zero mean here?
        return np.array(s, dtype=np.float32).reshape([self.L, self.L])

    def _append_record(self, record):
        with open(self.json_file, 'a') as f:
            json.dump(record, f)
            f.write(os.linesep)

    def dump_env_states(self):
        # get current timestamp
        total_steps = self.sim.get_total_steps()
        ep = self.sim.get_episode()
        # agent walk # steps in this episode
        ep_step_counters = self.sim.get_ep_step_counter()
        trajectory = self.sim.get_trajectory()
        if self.sim.get_accepted_length():
            loop_length = self.sim.get_accepted_length()[-1]
        else :
            loop_length = 0
        enclosed_area = self.calculate_area()
        update_times = self.sim.get_updated_counter()
        action_counters = self.sim.get_action_statistics()
        action_stats = [x / total_steps for x in action_counters]

        start_site = self.sim.get_start_point()
        acceptance = update_times * 100.0 / ep

        d = {
            'Episode': ep,
            'Steps'  : total_steps,
            'StartSite'  : start_site,
            'Trajectory': trajectory,
            'UpdateTimes': update_times,
            'AcceptanceRatio' : acceptance, 
            'LoopLength': loop_length,
            'EnclosedArea': enclosed_area,
            'ActionStats' : action_stats
        }
        
        self._append_record(d)


    ## ray test  (for: int, list)
    def conver_1Dto2D(self, input_1D):
        output_2D = None
        if type(input_1D) == int:
            output_2D = (int(input_1D/self.L), int(input_1D%self.L))
        elif type(input_1D) == list:
            output_2D = []
            for position in input_1D:
                output_2D.append((int(position/self.L), int(position%self.L)))
        return output_2D

    ## ray test
    def calculate_area(self):
        traj_2D = self.conver_1Dto2D(self.sim.get_trajectory())
        traj_2D_dict = {}
        for x, y in traj_2D:
            if x in traj_2D_dict:
                traj_2D_dict[x].append(y)
            else:
                traj_2D_dict[x] = [y]

        # check Max y_length
        y_position_list = []
        for y_list in traj_2D_dict.values():
            y_position_list = y_position_list + y_list
        y_position_list = list(set(y_position_list))
        max_y_length = len(y_position_list) -1

        area = 0
        for x in traj_2D_dict:
            diff = max(traj_2D_dict[x]) - min(traj_2D_dict[x])
            if diff > max_y_length:
                diff = max_y_length
            temp_area = diff - len(traj_2D_dict[x]) +1  ## avoid vertical straight line
            if temp_area > 0:
                area = area + temp_area

        return int(area)

    ## ray test
    def save_record_dict(self, length, area):
        def make_dict(target_dict, target_key):
            if target_key in target_dict:
                target_dict[target_key] = target_dict[target_key] +1
            else:
                target_dict[target_key] = 1
        make_dict(target_dict = self.record_dict[0], target_key = length)
        make_dict(target_dict = self.record_dict[1], target_key = area)

    ## ray test
    def hundred_test(self, accepted=False):
        hundred_size = 100
        if len(self.accepted_in_hundred_stack) < hundred_size:
            self.accepted_in_hundred_stack.append(accepted)
        elif accepted == True:
            self.accepted_in_hundred_stack[-1] = True
        elif len(self.accepted_in_hundred_stack) == hundred_size:
            self.accepted_in_hundred_stack = self.accepted_in_hundred_stack[1:] + [False]
        self.accepted_in_hundred = self.accepted_in_hundred_stack.count(True) 
