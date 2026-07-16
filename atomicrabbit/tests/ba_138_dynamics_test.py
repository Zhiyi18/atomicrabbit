# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:41:58 2026

@author: USER-01
"""

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
import qutip as qt
import numpy as np

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
    polarization = 0
    )

laser2 = lasers.Laser(
    frequency = transition_12.frequency, 
    intensity = 1e5,  
    linewidth = 4e6, 
    polarization = 0
    )

laser_list = [laser1, laser2]
state_list = [state1, state2, state3]

hyperfine_list = states.batch_hyperfine_split(state_list, I = 0, mu_I=0.1)
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
                  I = 0,
                  mu_I = 0.1,
                  B = 1e-4
    )

# driven_list = Ba137.driven_transition()
#for transition in driven_list:
    #print(f'Driven transition {transition.lower.spectroscopic_name} - {transition.upper.spectroscopic_name}')
    #print(f'Transition dipole moment = {transition.transition_dipole}')

N = 9

H, index = Ba137.build_hamiltonian()
H = qt.Qobj(H)

all_states = Ba137.add_zeeman_states()
ground_candidates = [s for s in all_states if s.parent.configuration == '6s2' and s.mF == 0]
print(ground_candidates)

psi0 = qt.basis(N, index[ground_candidates[0]])
rho0 = qt.ket2dm(psi0)

c_ops = []

t_end = 1e-8
n_steps = 500
tlist = np.linspace(0, t_end, num=n_steps)

result = qt.mesolve(H, 
                 rho0, 
                 tlist, 
                 c_ops, 
                 e_ops=[qt.basis(N, i) * qt.basis(N, i).dag() for i in range(N)])

#%%
print(len(ground_candidates))
for gc in ground_candidates:
    print(gc)
    print("F:", gc.F, "mF:", gc.mF, "config:", gc.parent.configuration)

idx = index[ground_candidates[0]]
print("psi0 index:", idx)
print("state at that index:", all_states[idx])
print("F of that state:", all_states[idx].F)

#%%
import matplotlib.pyplot as plt

groups = {}
for i, state in enumerate(all_states):
    F_val = getattr(state, 'F', 'unknown')
    groups.setdefault(F_val, []).append(i)

group_labels = [(indices, f"F={F}") for F, indices in groups.items()]

fig, axes = plt.subplots(len(group_labels), 1, figsize=(8, 3 * len(group_labels)), sharex=True)
if len(group_labels) == 1:
    axes = [axes]

colors = plt.cm.tab10.colors
for ax, (indices, title) in zip(axes, group_labels):
    for idx in indices:
        state = all_states[idx]
        F_val = getattr(state, 'F', '?')
        mF_val = getattr(state, 'mF', '?')
        ax.plot(tlist * 1e6, result.expect[idx], label = f"{state.parent.configuration} F={F_val}, mF={mF_val}", color=colors[idx % 10])
    ax.set_ylabel("Population")
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.set_ylim(-0.05, 1.05)

axes[-1].set_xlabel("Time (µs)")
plt.tight_layout()
plt.show()

#%%
# Group state indices by fine state (configuration + term symbol)
fine_groups = {}
for i, state in enumerate(all_states):
    fine_state = state.parent  # ZeemanState -> HyperfineState -> FineState
    key = fine_state.spectroscopic_name  # e.g. "6s2 1S_0"
    fine_groups.setdefault(key, []).append(i)

fig, axes = plt.subplots(len(fine_groups), 1, figsize=(9, 3.5 * len(fine_groups)), sharex=True)
if len(fine_groups) == 1:
    axes = [axes]

colors = plt.cm.tab20.colors
desired_order = ["6s2 1S_0", "6s6p 1P_1", "6s6d 1D_2"]
fine_groups_ordered = {key: fine_groups[key] for key in desired_order if key in fine_groups}

fig, axes = plt.subplots(len(fine_groups_ordered), 1, figsize=(9, 3.5 * len(fine_groups_ordered)), sharex=True)
if len(fine_groups_ordered) == 1:
    axes = [axes]

colors = plt.cm.tab20.colors

for ax, (label, indices) in zip(axes, fine_groups_ordered.items()):
    for j, idx in enumerate(indices):
        state = all_states[idx]
        F_val = state.F
        mF_val = state.mF
        ax.plot(tlist * 1e9, result.expect[idx],
                 label=f"F={F_val}, mF={mF_val}",
                 color=colors[j % 20])
    ax.set_ylabel("Population")
    ax.set_title(label)
    ax.legend(fontsize=7, ncol=2)
    ax.set_ylim(-0.05, 1.05)

axes[-1].set_xlabel("Time (ns)")
plt.tight_layout()
plt.show()