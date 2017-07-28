#from libicegame import SQIceGame, INFO
import gym
import gym_icegame
env = gym.make('IceGameEnv-v0')
import numpy as np

env.start(100)
for i in range(10):
    obs, r, done, info = env.step(np.random.randint(7))
