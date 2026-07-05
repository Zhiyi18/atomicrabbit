from dataclasses import dataclass
from typing import Optional
import numpy as np
from qutip import *
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
from atomicrabbit.physics import states

@dataclass(frozen=True)
class Transition:
    lower: states.FineState | states.HyperfineState | states.ZeemanState
    upper: states.FineState | states.HyperfineState | states.ZeemanState
    linewidth: Optional[float] = None
    lifetime: Optional[float] = None
    A: Optional[float] = None
    
    def __post_init__(self):
        if type(self.lower) is not type(self.upper):
            raise TypeError(
                f'Lower state type and upper state type mismatch when creating the transition.'
                f'Got {type(self.lower).__name__} and {type(self.upper).__name__}.'
            )
    
    
    @property
    def reduced_dipole(self):
        if self.A == None:
            raise TypeError(
                'Einstein coefficient A missing when creating the transition.'
                )
            
        d_ij = np.sqrt(
            3 * np.pi * epsilon_0 * hbar * c**3 * self.A
            / self.frequency()**3
            )
        return d_ij
    
    @property
    def hyperfine_factor(self):
        if isInstance(self.lower, FineState):
            return 1
        
        elif isInstance(self.lower, HyperfineState):
            hyperfine_factor = 
            return hyperfine_factor
        
        else isInstance(self.lower, ZeemanState):
            hyperfine_factor = 
            return hyperfine_factor
        
    @property
    def angular_factor(self):
        if isInstance(self.lower, FineState):
            angular_factor = 
            return angular_factor
        
        elif isInstance(self.lower, HyperfineState):
            angular_factor =
            return angular_factor
        
        else isInstance(self.lower, ZeemanState):
            angular_factor = 
            return angular_factor
        
    @property
    def transition_dipole(self):
        dipole = self. reduced_dipole
                * self. hyperfine_factor
                * self.angular_factor
        return dipole
        
        
    def frequency(self):
        return (self.upper.energy - self.lower.energy) * c *100
        
    def rabi_frequency(self, laser):
       E0 = laser.electric_field_amplitude()
       
       return self.dipole * E0 / hbar
    
    def B21(self):
        return self.A * c**3 / (8 * np.pi * h * self.frequency()** 3)
   
    def B12(self):
        g1 = self.lower.degeneracy()
        g2 = self.upper.degeneracy()
        return (g2 / g1) * self.B21()
    
    def detuning(self, laser):
        return self.frequency() - laser.frequency
    
    
        