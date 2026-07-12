from dataclasses import dataclass
from typing import Optional
import numpy as np
from qutip import *
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
from atomicrabbit.physics import states
from atomicrabbit.physics import lasers
from atomicrabbit.physics import transitions
import networkx as nx

@dataclass(frozen=True)
class Atom:
    species: str
    isotope: str
    transition_list: list[transitions.FineTransition]
    laser_list: list[lasers.Laser]
    I: Optional[float] = 0
    mu_I: Optional[float] = 0
    B: Optional[float] = 0
    
    # def load_transition(self, table):
        
    def add_fine_states(self):
        fine_states = set()
        for transition in self.transition_list:
            fine_states.add(transition.lower)
            fine_states.add(transition.upper)
        return list(fine_states)
    
    def add_zeeman_states(self):
        fine_states = self.add_fine_states()
        
        if self.I == 0:
            return states.batch_zeeman_split(fine_states, self.B)
        else:
            hyperfine_states = states.batch_hyperfine_split(fine_states, I=self.I, mu_I=self.mu_I)
            return states.batch_zeeman_split(hyperfine_states, self.B)
    
    # improve this later
    def zeeman_transitions(self):
        zeeman_transitions = []
        if self.I == 0:
            for transition in self.transition_list:
                lower_states = states.zeeman_split(transition.lower, self.B)
                upper_states = states.zeeman_split(transition.upper, self.B)
                for state1 in lower_states:
                    for state2 in upper_states:
                        zeeman_transitions.append(transitions.Transition(
                            lower = state1,
                            upper = state2,
                            A = transition.A,
                            linewidth = transition.linewidth,
                            lifetime = transition.lifetime
                            ))
                        
        else:
            for transition in self.transition_list:
                lower_states = states.batch_zeeman_split(states.hyperfine_split(transition.lower, I = self.I, mu_I = self.mu_I), self.B)
                upper_states = states.batch_zeeman_split(states.hyperfine_split(transition.upper, I = self.I, mu_I = self.mu_I), self.B)
                for state1 in lower_states:
                    for state2 in upper_states:
                        zeeman_transitions.append(transitions.Transition(
                            lower = state1,
                            upper = state2,
                            A = transition.A,
                            linewidth = transition.linewidth,
                            lifetime = transition.lifetime
                            ))
                        
        return zeeman_transitions
    
    def driven_transitions(self):
        # add situations where the transitions are not zeeman states later
        driven_transition = []
        zeeman_transitions = self.zeeman_transitions()
        for zeeman_transition in zeeman_transitions:
            for laser in self.laser_list:
                dipole = zeeman_transition.transition_dipole(laser)
                if dipole != 0:
                    driven_transition.append(transitions.DrivenTransition(
                        transition = zeeman_transition,
                        laser = laser,
                        dipole = dipole
                        ))
                    print(f'Driven transition {zeeman_transition.lower.spectroscopic_name} - {zeeman_transition.upper.spectroscopic_name}')
                    print(f'Transition dipole = {dipole}')
                    print(f'Laser = {laser}')                
                        
        return driven_transition
    
    @property
    def system_size(self):
        # add non-zeeman states later too
        # this thing calls zeeman_transitions twice and might slow the calculations down!
        return len(self.zeeman_transitions())
    
    def atom_graph(self):
        G = nx.Graph()
        G.add_nodes_from(self.add_zeeman_states())
        
        for driven in self.driven_transitions():
            lower = driven.transition.lower
            upper = driven.transition.upper
            Omega = driven.rabi_frequency
            detuning = driven.detuning
            G.add_edge(lower, upper, Omega=Omega, detuning=detuning)
        
        return G
    
    def assign_rotating_frame_energies(self):
        # I don't know how this works and will come back later
        G = self.atom_graph()
        energies = {}
        for component in nx.connected_components(G):
            ref = min(component, key=lambda s: s.energy_omega)
            sub_energies = {ref: ref.energy_omega}
            
            for u, v in nx.bfs_edges(G, ref):
                detuning = G[u][v]['detuning']
                sub_energies[v] = sub_energies[u] - hbar * detuning
                energies.update(sub_energies)
    
        return energies
    
    def build_hamiltonian(self):
        states = self.add_zeeman_states()
        n = len(states)
        index = {s: i for i, s in enumerate(states)}
        H = np.zeros((n, n), dtype=complex)
        
        energies = self.assign_rotating_frame_energies()
        for state, E in energies.items():
            H[index[state], index[state]] = E
        
        for driven in self.driven_transitions():
            i, j = index[driven.transition.lower], index[driven.transition.upper]
            H[i, j] = driven.rabi_frequency / 2
            print(f'rabi frequency = {driven.rabi_frequency}')
            H[j, i] = np.conj(H[i, j])
        
        return H, index
                
                    

