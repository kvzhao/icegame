import logging
from gym.envs.registration import register

logger = logging.getLogger(__name__)

register (
    id = 'SAWCanvasEnv-v0',
    entry_point='gym_saw.envs:SAWCanvasEnv',
    kwargs={
        'L' : 32,
    },
    )