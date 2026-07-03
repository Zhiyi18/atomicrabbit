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
    
    
    @property
    #def dipole_moment(self):
        
        
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
    
        