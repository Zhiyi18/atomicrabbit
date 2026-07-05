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
    A_hyperfine = -109.5,
    B_hyperfine = 0
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

'''
state3_split = states.hyperfine_split(state3, I = 3/2)
print(state3_split)
for state in state3_split:
    print(state.energy_cm)
    
state3_zeeman_split = states.zeeman_split(state3, B=1e-4, gF=2)
for state in state3_zeeman_split:
    print(state.energy_cm)
'''

state_list = [state1, state2, state3]
hyperfine_list = states.batch_hyperfine_split(state_list, I = 3/2, mu_I=0.1)
for state in hyperfine_list:
    print(state.energy_cm)
    
zeeman_list = states.batch_zeeman_split(state_list, B=1e-4)
for state in zeeman_list:
    print(state.energy_cm)

full_splitting = states.batch_zeeman_split(hyperfine_list, B=1e-4)
for state in full_splitting:
    print(state.energy_cm)
print(len(full_splitting))


