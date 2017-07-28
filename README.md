# icegame

## intro
Designed as an environments interacting with spin ice system.

## Compile & Install 

### compile

```
sh compile.sh
```

will generate icegame.so in build/src. please move icegame.so to where you executing experiment

```
mv build/src/icegame.so /where/you/run/icegame
```

### gym-icegame

Install python interface 

```
python setup.py install
```

which depends on openai gym.


## Game Scenario
Draw a proposed loop, then summit.

### actions
* Directional action: (up, down, left, right, next up, next down) 6 operations in total

Two options: 
* Metropolis button: 
* Loop Auto-detection: 

function called by env.step(action_index)
* step()

### rewards

* Step-wise reward
```
    r = scale * (icemove_w * rets[0] + energy_w * rets[1] + defect_w * rets[2] + baseline)
```
* Accepted reward
```
    r = +1.0 * (loop_length / 4.0)
```

```
rets[0]: +1/-1 accpet/reject
rets[1]: difference of totoal mean energy 
rets[2]: difference of defect density
rets[3]: ratio of configuration difference
```

### observations
stacked scene
* spin configuration
* trajectory
* energy map
* defect map
in format of `NHWC`


## Callable interface from libicegame
List in sqice.hpp

## TODO
* Long and short loop detection mechanism
* <s>Reward design</s>
* Area reward
* <s>Save loop config when accepted</s>
* <s>Fix ambiguous `rets`</s>
* demo codes
* Analysis tools

