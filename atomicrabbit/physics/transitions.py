from dataclasses import dataclass
from typing import Optional
import numpy as np
from qutip import *
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
from sympy.physics.wigner import wigner_3j, wigner_6j
from sympy import N
from atomicrabbit.physics import states
from atomicrabbit.physics import lasers

@dataclass(frozen=True)
class FineTransition:
    lower: states.FineState
    upper: states.FineState
    A: float
    linewidth: float | None = None
    lifetime: float | None = None
    
@dataclass(frozen=True)
class Transition(FineTransition):
    lower: states.FineState | states.HyperfineState | states.ZeemanState
    upper: states.FineState | states.HyperfineState | states.ZeemanState

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
                'Einstein coefficient is required to calculate the transition dipole element.'
                )
            
        d_ij = np.sqrt(
            3 * np.pi * epsilon_0 * hbar * c**3 * self.A
            / self.frequency**3
            )
        return d_ij
    
    @property
    def hyperfine_factor(self):
        if isinstance(self.lower, states.FineState):
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
                    * np.sqrt((2 * Fp + 1) * (2 * F + 1))
                    * hf
                )
        

    def angular_factor(self, laser):
        q = laser.polarization
        if q not in(-1, 0, 1):
            raise ValueError(
                'Laser polarisation must be -1, 0, or 1'
                )
            
        elif isinstance(self.lower, states.FineState):
            Fp = self.upper.J
            mp = self.upper.mJ
            F = self.lower.J
            m = self.lower.mJ
            exponent = Jp - mp
        
        elif isinstance(self.lower, states.HyperfineState):
            raise ValueError(
                'HyperfineState has no mF. Use ZeemanState to calculate angular factors.'
                )
        
        elif isinstance(self.lower, states.ZeemanState):
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
        return (self.upper.energy_cm - self.lower.energy_cm) * c *100
        
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
        return self.frequency - laser.frequency
    
@dataclass(frozen=True)
class DrivenTransition:
    transition: Transition
    laser: lasers.Laser
    dipole: float
    
    @property
    def lower(self):
        return self.transition.lower
    
    @property
    def upper(self):
        return self.transition.upper
    
    @property
    def rabi_frequency(self):
        return self.transition.rabi_frequency(self.laser)
    
    @property
    def detuning(self):
        return self.transition.detuning(self.laser)
    

