from dataclasses import dataclass
from typing import Optional
import numpy as np
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants

'''
example use:
    
d52 = FineState(...)

F3 = HyperfineState(
    parent=d52,
    F=3
)

mF2 = ZeemanState(
    parent=F3,
    m=2
)

'''
mu_B = physical_constants["Bohr magneton"][0]

@dataclass(frozen=True)
class FineState:
    configuration: str
    term: str
    L: float
    S: float
    J: float
    energy_cm: float
    parity: str
    I: Optional[float] = None
    gJ: Optional[float] = None
    
    @property
    def g(self):
       return self.gJ
   
    @property
    def degeneracy(self):
        return 2 * self.J +1
    
    @property
    def term_symbol(self):
        L_letters = {
            0: "S",
            1: "P",
            2: "D",
            3: "F",
            4: "G",
        }

        multiplicity = int(2 * self.S + 1)
        L_letter = L_letters.get(self.L, "?")

        return f"{multiplicity}{L_letter}_{self.J}"
    
    @property
    def spectroscopic_name(self):
        return f"{self.configuration} {self.term_symbol()}"


@dataclass(frozen=True)
class HyperfineState:
    parent: FineState
    F: float
    gF: float
    A: float
    B: float = 0

    @property
    def J(self):
        return self.fine_state.J

    @property
    def I(self):
        return self.fine_state.I
    
    @property
    def g(self):
        return self.gF
    
    @property
    def hyperfine_shift(self):
        '''
        calculate hyperfine shift
        '''
        
        J = self.J
        I = self.I
        F = self.F
        A = self.A
        B = self.B

        C = F * (F + 1) - I * (I + 1) - J * (J + 1)
    
        # Magnetic dipole term
        delta_E = 0.5 * A * C
    
        # Electric quadrupole term
        if B != 0 and I > 0.5 and J > 0.5:
            numerator = (
                1.5 * C * (C + 1)
                - 2 * I * (I + 1) * J * (J + 1)
            )
    
            denominator = (
                I * (2 * I - 1)
                * J * (2 * J - 1)
            )
    
            delta_E += (B / 4) * numerator / denominator
    
        return delta_E
    
    @property
    def energy_cm(self):
        return (
            self.parent.energy_cm
            + self.hyperfine_shift()
        )

    
@dataclass(frozen=True)
class ZeemanState:
    parent: FineState | HyperfineState
    mF: float
    B : float   # Magnetic field strength

    @property
    def zeeman_shift(self):
        delta_E = mu_B * self.parent.g * self.mF * self.B
        return delta_E
    
    @property
    def energy_cm(self):
        return (
            self.parent.energy_cm
            + self.zeeman_shift()
        )
    