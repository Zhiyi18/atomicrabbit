from dataclasses import dataclass
from typing import Optional
import numpy as np
from qutip import *
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
from sympy.physics.wigner import wigner_3j, wigner_6j
from sympy import N
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
        if self.A is None:
            raise ValueError(
                'Einstein coefficient A missing when creating the transition.'
                )
            
        d_ij = np.sqrt(
            3 * np.pi * epsilon_0 * hbar * c**3 * self.A
            / self.frequency()**3
            )
        return d_ij
    
    @property
    def hyperfine_factor(self):
        if isinstance(self.lower, FineState):
            return 1
        
        else:
            Jp = self.upper.J
            Fp = self.upper.F
            I = self.lower.I
            J = self.lower.J
            F = self.lower.F
            
            hf = float(N(wigner_6j(Jp, Fp, I, F, J, 1)))
            
            return (
                    (-1) ** (Jp + I + F + 1)
                    * sqrt((2 * Fp + 1) * (2 * F + 1))
                    * hf
                )
        

    def angular_factor(self, laser):
        q = laser.polarisation
        if laser.polarisation not in(-1, 0, 1):
            raise ValueError(
                'Laser polarisation must be -1, 0, or 1'
                )
            
        elif isinstance(self.lower, FineState):
            Fp = self.upper.J
            mp = self.upper.mJ
            F = self.lower.J
            m = self.lower.mJ
            exponent = Jp - mp
        
        elif isinstance(self.lower, HyperfineState):
            raise ValueError(
                'HyperfineState has no mF. Use ZeemanState to calculate angular factors.'
                )
        
        elif isinstance(self.lower, ZeemanState):
            Fp = self.upper.F
            mp = self.upper.mF
            F = self.lower.F
            m = self.lower.mF
            exponent = Fp - mp
            
        else:
            raise TypeError(
                f'Unsupported state type: {type(self.lower).__name__}'
                )
            
        cg = float(N(wigner_3j(Fp, 1, F, -mp, q, m)))
        return (
                (-1) ** exponent
                * cg
        )
        

    def transition_dipole(self, laser):
        dipole = (
            self. reduced_dipole
            * self. hyperfine_factor
            * self.angular_factor(laser)
        )
        return dipole
        
    @property
    def frequency(self):
        return (self.upper.energy - self.lower.energy) * c *100
        
    def rabi_frequency(self, laser):
       E0 = laser.electric_field_amplitude()
       
       return self.transition_dipole(laser) * E0 / hbar
    
    def B21(self):
        return self.A * c**3 / (8 * np.pi * h * self.frequency()** 3)
   
    def B12(self):
        g1 = self.lower.degeneracy()
        g2 = self.upper.degeneracy()
        return (g2 / g1) * self.B21()
    
    def detuning(self, laser):
        return self.frequency() - laser.frequency
    
    
        