import matplotlib.pyplot as plt
import sys, os
import json
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Analyise history from Icegame environment.')

parser.add_argument('--logdir', type=str, default='env_history.json',
                   help='Json log file path.')

args = parser.parse_args()

def read_json(filename):
    with open(filename, 'r') as json_file:
        d_list = [json.loads(line) for line in json_file]
        print ('Log contains {} timestamps.'.format(len(d_list)))
    return d_list

def site2coord(sites, L):
    coords = None
    if type(sites) == int:
        coords = (int(sites/L), int(sites%L))
        # return tuple
    elif type(sites) == list:
        coords = [] # return list
        for s in sites:
            coords.append((int(s/L), int(s%L)))
    return coords 

def plot_loops(loops2d):
    x, y = zip(*loops2d)
    plt.scatter(x, y, vmax=32)
    plt.title('Loop configurations')
    plt.savefig('loop_config.png')
    plt.show()

# then we should carefully parse list

data = read_json(args.logdir)
# extract data information
'''
    Several items stored in data
    * Episode
    * Steps
    * Trajectory (list)
    * EnclosedArea
    * ActionStats (list)
    * LoopLength
    * UpdateTimes
    * AcceptanceRatio
    * StartSite
'''

episodes = []
steps = []
updates = []
lengths = []
areas = []
accept_ratio = []
action_stats = []
loops = []
loops2d = []

indices = list(range(len(data)))

## extract data
for d in data:
    episodes.append(d['Episode'])
    steps.append(d['Steps'])
    lengths.append(d['LoopLength'])
    areas.append(d['EnclosedArea'])
    accept_ratio.append(d['AcceptanceRatio'])
    action_stats.append(d['ActionStats'])
    loops.append(d['Trajectory'])
    loops2d.extend(site2coord(d['Trajectory'], 32))

plot_loops(loops2d)

plt.plot(indices, accept_ratio)
plt.title('Acceptance Ratio')
plt.ylabel('updated times / total episodes')
plt.savefig('acc_ratio.png')
plt.show()

plt.plot(indices, lengths, linewidth=1.5)
plt.plot(indices, areas, linewidth=1.5)
axes = plt.gca()
axes.set_ylim([-1, 32])
plt.legend(['Accepted Length', 'Enclosed Area'])
plt.title('Accepted Loop Topology')
plt.savefig('loop_topo.png')
plt.show()