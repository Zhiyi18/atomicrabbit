from dataclasses import dataclass
from typing import Optional
from typing import Literal
import numpy as np
from scipy.constants import h, hbar, epsilon_0, c, k, physical_constants
from atomicrabbit.tools import function_generator
    
    
@dataclass(frozen=True, kw_only=True)
class Laser:
    """
    Properties
    ---------
    frequency: Hz
    intensity: W * m^-2
    linewidth: Hz
    ---------
    All axes are given in Cartesian coordinates.
    If not specified, the default propagation direction is in z axis direction.
    -----------------------------
    Polarization is given in Spherical basis or Jones basis.
    Spherical basis: +/- 1 for sigma +/- polarization, and 0 for pi polarization.
    Jones basis: Polarization represented in Jones vector.
    -----------------------------
    

    """
    frequency: float
    intensity: float
    linewidth: float
    propagation_direction: Optional[tuple[float, float, float]] = [0, 0, 1]
    polarization_basis: Literal["spherical", "jones"]
    polarization: int | tuple[complex, ...]
    quantization_axis: Optional[tuple[float, float, float]] = None

    """
    Add pulse shape later
    There should be pulse objects like pi pulse, gaussian pulse, etc...
    Class Pulse is to be added in function_generator(maybe move it to physics?)
    pulse_shape: Optional[function_generator.Pulse]
    """
    
    def __post_init__(self):
            if self.polarization_basis == "spherical":
                if self.polarization not in (-1, 0, 1):
                    if not isinstance(self.polarization, tuple) or len(self.polarization) != 3:
                        raise ValueError(
                            "Spherical polarization must be q=-1/0/+1 or a 3-component tuple."
                        )
                    
                if self.quantization_axis is not None:
                    print(
                        "The quantization axis information is not used because the polarization is given in spherical basis"
                    )
    
            elif self.polarization_basis == "jones":
                if not isinstance(self.polarization, tuple) or len(self.polarization) != 2:
                    raise ValueError("Jones polarization must be a 2-component tuple.")
                
                if self.quantization_axis is None:
                    raise ValueError("Jones polarization requires a quantization_axis.")
    
            else:
                raise ValueError("polarization_basis must be 'spherical' or 'jones'.")
            
    def electric_field_amplitude(self):
        return np.sqrt(2 * self.intensity / (c * epsilon_0))

    # Generating a list of pulse amplitudes with time list
    def pulse(self, pulse_shape, t):
        if pulse_shape is 'cw':
            return self.electric_field_amplitude * t

    def spherical_components(self, quantization_axis=None):
        if self.polarization_basis == "spherical":
            # Convert q shorthand or normalize the supplied 3-vector.
            
        else:
            axis = quantization_axis or self.quantization_axis
            # Convert Jones components using propagation direction and axis.
            

        


    
    