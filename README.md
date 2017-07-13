# icegame

## intro
Designed as an environments interacting with spin ice system.

## Compile & Install 

### src

```
make -j4; make icegame
```

will generate libicegame.so which should be put under your execution folder.

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

function called by env.step(action_index), here we provide severl step()
* step()
* step_auto()
* step_binary()

### rewards

* Step-wise reward
```
    r = - 10.0 * (energy difference + defect density)
```
* Accepted reward
```
    r = +1.0
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
* Short loop flip trajectory
* Fix ambiguous `rets`
* demo codes

