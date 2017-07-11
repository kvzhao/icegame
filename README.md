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


## TODO
* Long and short loop detection mechanism
* <s>Reward design</s>
* Area reward
* Save loop config when accepted
* Short loop flip trajectory

