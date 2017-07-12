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


### Game Scenario

#### actions
* Directional action (up, down, left, right, next up, next down) 6 operations

Two options: 
* Metropolis button
* Auto-detection

#### rewards

* Step-wise reward
```
    r = - 10.0 * (energy difference + defect density)
```
* Accepted reward
```
    r = +1.0
```

## TODO
* Long and short loop detection mechanism
* <s>Reward design</s>
* Area reward
* Save loop config when accepted
* Short loop flip trajectory
* <s>Save loop configuration</s>

