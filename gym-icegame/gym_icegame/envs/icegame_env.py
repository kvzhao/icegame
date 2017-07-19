from __future__ import division
import gym
from gym import error, spaces, utils, core

#from six import StingIO
import sys
import numpy as np
import random
from libicegame import SQIceGame, INFO

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
                                  7 :   'noop',
                                  })

        self.index_mapping = dict({
                                  'right': 0,
                                  'down' : 1,
                                  'left' : 2,
                                  'up' : 3,
                                  'lower_next' : 4,
                                  'upper_next' : 5,
                                  'metropolis' : 6,
                                  'noop'       : 7,
                                  })


        ### action space and state space
        self.action_space = spaces.Discrete(len(self.name_mapping))
        self.observation_space = spaces.Box(low=-1.0, high=1.0, shape=(self.L, self.L, 4))
        self.reward_range = (-1, 1)

        # output file
        self.ofilename = 'loop_sites.log'
        # render file
        self.rfilename = 'loop_renders.log'
        self.stacked_axis = 2

        ## counts reset()
        self.episode_counter = 0

    def step(self, action):
        terminate = False
        reward = 0.0
        obs = None
        rets = [0.0, 0.0, 0.0, 0.0]
        metropolis_executed = False

        if (action == 6):
            self.sim.flip_trajectory()
            rets = self.sim.metropolis()
            metropolis_executed = True
        elif (0 <= action < 6) :
            rets = self.sim.draw(action)

        if (metropolis_executed):
            if rets[0] > 0 and rets[3] > 0:
                print ('ACCEPTS!')
                self.sim.update_config()
                loop_length = self.sim.get_accepted_length()[-1]
                reward = 1.0 * (loop_length / 4.0) # reward with different length by normalizing with len 4 elements
                with open(self.ofilename, 'a') as f:
                    f.write('{}\n'.format(self.sim.get_trajectory()))
                    print ('\tSave loop configuration to file')
                print ('\tTotal accepted number = {}'.format(self.sim.get_updated_counter()))
                print ('\tAccepted loop len = {}'.format(loop_length))
                print ('\tTotal accepted loop length = {}'.format(self.sim.get_accepted_length()))
                self.render()
                self.sim.clear_buffer()
            else:
                self.sim.clear_buffer()
                # Avoid running metropolis at start
                if (rets[3] == 0.0):
                    reward = -0.8
            # reset or update
        else:
            reward = self._stepwise_weighted_returns(rets)
            # as usual

        obs = self.get_obs()
        ## add timeout mechanism?

        return obs, reward, terminate, rets

    ## DEPRECATED!! 
    def step_auto(self, action):
        terminate = False
        reward = 0.0
        obs = None
        rets = [0.0, 1.0, 1.0]
        metropolis_executed = False

        '''
            1. detection stage
                if dd & de == 0:
                    flip long loop
                    run metropolis
                if short loop detected:
                    flip short loop
                    run metropolis
                else:
                    walk
            2. reward eval stage
                if metropolis executed:
                    reward with accpetance or not
                if as usual:
                    reward is calculated by dE and dD
        '''
        
        # restrict action from 0 to 4: {up, down, left, right}
        if (0<= action < 4):
            # ignore other action idxes
            rets = self.sim.draw(action)
        else:
            print ('Accept only directional action: up/down/left/right!')

        dE = rets[1]
        dD = rets[2]
            
        ### Think about
        if (dE == 0.0 and dD == 0.0):
            self.sim.flip_trajectory()
            rets = self.sim.metropolis()
            metropolis_executed = True
        #TODO: Short loop detection
        
        # EVALUATION
        if (metropolis_executed):
            if rets[0] > 0 and rets[3] > 0:
                print ('ACCEPTS!')
                reward = 1.0
                self.sim.update_config()
                print (self.sim.get_trajectory())
                with open(self.ofilename, 'a') as f:
                    f.write('{}\n'.format(self.sim.get_trajectory()))
                    print ('\tSave loop configuration to file')
                print ('\tTotal accepted number = {}'.format(self.sim.get_updated_counter()))
                print ('\tAccepted loop length = {}'.format(self.sim.get_accepted_length()))
                self.sim.clear_buffer()
            else:
                self.sim.clear_buffer()
                if (rets[3] ==  0):
                    reward = -0.8
            # reset or update
        else:
            reward = self._stepwise_weighted_returns(rets)
            # as usual
        # RETURN
        ## check timeout!
        if (self.timeout()):
            terminate = True

        obs = self.get_obs()
        return obs, reward, terminate, rets

    # Start function used for agent learing
    def start(self, init_site):
        init_agent_site = self.sim.start(init_site)
        assert(init_site == init_agent_site)

    def reset(self):
        ## clear buffer and set new start of agent
        self.sim.clear_buffer()
        self.start(rnum(self.N))
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

    def _stepwise_weighted_returns(self, rets):
        icemove_w = 0.000
        energy_w = -1.0
        defect_w = 0.0
        baseline = 0.009765625 ## 1 / 1024
        scaling = 2.0
        return (icemove_w * rets[0] + energy_w * rets[1] + defect_w * rets[2] + baseline) * scaling

    def sample_icemove_action_index(self):
        return self.sim.icemove_index()

    def render(self, mapname ='traj', mode='ansi', close=False):
        #of = StringIO() if mode == 'ansi' else sys.stdout
        #print ('Energy: {}, Defect: {}'.format(self.sqice.cal_energy_diff(), self.sqice.cal_defect_density()))
        s = None
        if (mapname == 'traj'):
            s = self._transf2d(self.sim.get_canvas_map())
        elif (mapname == 'state'):
            s = self._transf2d(self.sim.get_state_t_map())
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
            screen += '|\n'
        screen += '\t+' + self.L * '---' + '+\n'
        #sys.stdout.write(screen)
        with open(self.rfilename, 'a') as f:
            f.write('Episode: {}\n'.format(self.episode_counter))
            f.write('{}\n'.format(screen))

    def get_obs(self):
        config_map = self._transf2d(self.sim.get_state_t_map())
        canvas_map = self._transf2d(self.sim.get_canvas_map())
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
