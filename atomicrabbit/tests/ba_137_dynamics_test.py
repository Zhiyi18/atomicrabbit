# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 19:38:36 2026

@author: USER-01
"""

# Add naming convention rules!
# Class: AtomicRabbit
# Instance: atomicRabbit
# Method: atomic_rabbit
from atomicrabbit.physics import states

state1 = states.FineState(
    configuration = '6s2',
    L = 0,
    S = 0,
    J = 0,
    energy_cm = 0,
    parity = '',
    )

state2 = states.FineState(
    configuration = '6s6p',
    L = 1,
    S = 0,
    J = 1,
    energy_cm = 18060.261,
    parity = '',
    A_hyperfine = -109.5,
    B_hyperfine = 0
    )

state3 = states.FineState(
    configuration = '6s6d',
    L = 2,
    S = 0,
    J = 2,
    energy_cm = 30236.826,
    parity = '',
    A_hyperfine = -109.5,
    B_hyperfine = 0
    )

state3_split = states.hyperfine_split(state3, I = 3/2)
print(state3_split)
for state in state3_split:
    print(state.energy_cm)



