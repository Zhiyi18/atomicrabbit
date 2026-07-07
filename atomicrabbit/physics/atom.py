from dataclasses import dataclass
from typing import Optional
import numpy as np
from qutip import *
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
from atomicrabbit.physics import states
from atomicrabbit.physics import lasers
from atomicrabbit.physics import transitions

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
        return states.batch_zeeman_split(fine_states)
    
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
                            

    def driven_transition(self):
        driven_transition = []
        zeeman_transitions = self.zeeman_transitions()
        for zeeman_transition in zeeman_transitions:
            for laser in self.laser_list:
                dipole = zeeman_transition.transition_dipole(laser)
                if dipole != 0:
                    driven_transition.append(zeeman_transition)
                    print(f'Driven transition {zeeman_transition.lower.spectroscopic_name} - {zeeman_transition.upper.spectroscopic_name}')
                    print(f'Transition dipole = {dipole}')
                        
        return driven_transition