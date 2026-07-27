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
        
    def basis(self):
        states = self.add_zeeman_states()
        index = {s: i for i, s in enumerate(states)}
        return states, index
    
    '''
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
    '''
    def magnetic_sublevels(self, state):
        if self.I == 0:
            return states.zeeman_split(state, self.B)
    
        return states.batch_zeeman_split(
            states.hyperfine_split(
                state,
                I=self.I,
                mu_I=self.mu_I,
            ),
            self.B,
        )

    def zeeman_transitions(self):
        zeeman_transitions = []
    
        for transition in self.transition_list:
    
            lower_states = self.magnetic_sublevels(transition.lower)
            upper_states = self.magnetic_sublevels(transition.upper)
    
            for lower in lower_states:
                for upper in upper_states:
                    zeeman_transitions.append(
                        transitions.Transition(
                            lower=lower,
                            upper=upper,
                            A=transition.A,
                            linewidth=transition.linewidth,
                            lifetime=transition.lifetime,
                        )
                    )

        return zeeman_transitions
    
    def branching_ratios(self):
        """Returns normalized branching ratios for every Zeeman transition."""
    
        ratios = {}
    
        for transition in self.transition_list:
    
            lower_states = self.magnetic_sublevels(transition.lower)
            upper_states = self.magnetic_sublevels(transition.upper)
    
            # Collect all allowed decay amplitudes
            amplitudes = {}
    
            for upper in upper_states:
                total = 0.0
    
                for lower in lower_states:
    
                    decay = transitions.Transition(
                        lower=lower,
                        upper=upper,
                        A=transition.A,
                        linewidth=transition.linewidth,
                        lifetime=transition.lifetime,
                    )
    
                    amp2 = abs(decay.branching_amplitude)**2
    
                    if amp2 > 0:
                        amplitudes[decay] = amp2
                        total += amp2
    
                # Normalize excited state's decay channels
                for decay, amp2 in amplitudes.items():
                    if decay.upper == upper:
                        ratios[decay] = amp2 / total
    
        return ratios

    def driven_transitions(self):
        # Add situations where the transitions are not zeeman states later
        driven_transition = []
        zeeman_transitions = self.zeeman_transitions()
        for zeeman_transition in zeeman_transitions:
            for laser in self.laser_list:
                dipole = zeeman_transition.transition_dipole(laser)
                detuning = zeeman_transition.detuning(laser)
                
                # Dropping the dipole-forbidden terms
                if dipole != 0:
                    
                    # Dropping the far-off-resonant transitions
                    if abs(detuning) <= 1e13:
                        driven_transition.append(transitions.DrivenTransition(
                            transition = zeeman_transition,
                            laser = laser,
                            dipole = dipole
                            ))
                        
                    print(f'Driven transition {zeeman_transition.lower.spectroscopic_name} - {zeeman_transition.upper.spectroscopic_name}')
                    print(f'Transition dipole = {dipole}')
                    print(f'Laser = {laser}')                
                        
        return driven_transition
    
    def driven_transition_debug(self):
        for driven in self.driven_transitions():
            print(
                f"{driven.lower.spectroscopic_name}"
                f" --> "
                f"{driven.upper.spectroscopic_name}"
                f"   q={driven.laser.polarization}"
                f"   Ω={driven.rabi_frequency:.3e}"
                f"   Δ={driven.detuning:.3e}"
            )
        
    @property
    def system_size(self):
        # add non-zeeman states later too
        # this thing calls zeeman_transitions twice and might slow the calculations down!
        return len(self.zeeman_transitions())
    
    def atom_graph(self):
        G = nx.Graph()
        G.add_nodes_from(self.add_zeeman_states())
                    
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
                sub_energies[v] = sub_energies[u] - 2 * np.pi * detuning
                # sub_energies[v] = sub_energies[u] - hbar * detuning
            energies.update(sub_energies)
    
        return energies
    
    def build_hamiltonian(self):
        # Now this part is shit because the diagonal terms and
        # the off-diagonal terms are added separately from the 
        # transitions and the atom graph
        
        states, index = self.basis()
        n = len(states)
        H = np.zeros((n, n), dtype=complex)
        
        # Adding diagonal terms
        energies = self.assign_rotating_frame_energies()
        for state, E in energies.items():
            # Dropping the off-resonant optical detunings
            if abs(E) <= 1e10:
                H[index[state], index[state]] = E
                print(f'diagonal detuning = {E} added')
        
        # Adding off-diagonal terms
        for driven in self.driven_transitions():
            i, j = index[driven.transition.lower], index[driven.transition.upper]
            H[i, j] += driven.rabi_frequency / 2
            print(f'rabi frequency = {driven.rabi_frequency}')
            H[j, i] = np.conj(H[i, j])

        return H, index
                
    def build_collapse_operators(self):

        states, index = self.basis()
    
        branching = self.branching_ratios()
    
        c_ops = []
    
        for transition, ratio in branching.items():
    
            gamma = transition.A * ratio
    
            if gamma == 0:
                continue
    
            lower = index[transition.lower]
            upper = index[transition.upper]
    
            c = (
                np.sqrt(gamma)
                * basis(len(states), lower)
                * basis(len(states), upper).dag()
            )
    
            c_ops.append(c)
    
        return c_ops
    
