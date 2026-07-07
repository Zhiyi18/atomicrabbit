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
from atomicrabbit.physics import lasers
from atomicrabbit.physics import transitions
from atomicrabbit.physics import atom

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

transition_12 = transitions.Transition(state1, state2, linewidth=20e6, A=1.19e8)
transition_23 = transitions.Transition(state2, state3, linewidth=20e6, A=5.00e6)
transition_13 = transitions.Transition(state1, state3, linewidth=20e6, A=1e5)
transition_list = [transition_12, transition_23, transition_13]

laser1 = lasers.Laser(
    frequency = transition_12.frequency, 
    intensity = 1e5,  
    linewidth = 4e6, 
    polarization = 1
    )

laser2 = lasers.Laser(
    frequency = transition_12.frequency, 
    intensity = 1e5,  
    linewidth = 4e6, 
    polarization = 1
    )

laser_list = [laser1, laser2]

'''
laser3 = lasers.Laser(
    frequency = c / 405e-9,
    intensity = 1e10, 
    linewidth=4e6
    )
'''

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
#for state in hyperfine_list:
    #print(state.energy_cm)
    
zeeman_list = states.batch_zeeman_split(state_list, B=1e-4)
#for state in zeeman_list:
    #print(state.energy_cm)

full_splitting = states.batch_zeeman_split(hyperfine_list, B=1e-4)
for state in full_splitting:
    print(state.energy_cm)
print(len(full_splitting))


Ba137 = atom.Atom(species = 'Ba I',
                  isotope = '137',
                  transition_list = transition_list,
                  laser_list = laser_list,
                  I = 3/2,
                  mu_I = 0.1,
                  B = 1e-4
    )

driven_list = Ba137.driven_transition()
#for transition in driven_list:
    #print(f'Driven transition {transition.lower.spectroscopic_name} - {transition.upper.spectroscopic_name}')
    #print(f'Transition dipole moment = {transition.transition_dipole}')
