from dataclasses import dataclass
from typing import Optional
import numpy as np
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
    
    
@dataclass(frozen=True)
class Laser:
    frequency: float
    intensity: float
    linewidth: float
    polarization: Optional[int] = None
    
    def __post_init__(self):
        if self.polarization not in (-1, 0, 1):
            raise ValueError(
                "Laser polarization must be -1, 0 or +1."
            )
            
    def electric_field_amplitude(self):
        return np.sqrt(2 * self.intensity / (c * epsilon_0))
    
    