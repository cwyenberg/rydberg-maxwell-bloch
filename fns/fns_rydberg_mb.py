import numpy as np
import scipy as sp
import qutip as qtp
import re
import pandas as pd

def parse_spec_state(state_str):
    """
    Parses a spectroscopic state string and returns the quantum numbers.
    Args:
        state_str (str): A string representing the spectroscopic state in the format "nlj",
                         where 'n' is the principal quantum number (integer),
                         'l' is the orbital type (one of 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k'),
                         and 'j' is the total angular momentum quantum number (float).
    Returns:
        tuple: A tuple containing three elements:
               - n (int): The principal quantum number.
               - l (int): The orbital quantum number.
               - j (float): The total angular momentum quantum number.
    Raises:
        ValueError: If the input string does not match the expected format.
    """
    
    # Define a dictionary to convert orbital symbols to quantum numbers
    l_mapping = {
        's': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6, 'k': 7
    }
    
    # Regular expression to extract n, l, and j
    match = re.match(r"(\d+)([spdfghik])([\d.]+)", state_str)
    if not match:
        raise ValueError(f"Invalid spectroscopic state string: {state_str}")
    
    # Extract matched groups
    n = int(match.group(1))  # Principal quantum number
    l_symbol = match.group(2)  # Orbital type
    j = float(match.group(3))  # Total angular momentum quantum number
    
    # Convert orbital symbol to quantum number l
    l = l_mapping[l_symbol]
    
    return n, l, j


def build_state_str(n: int, l: int, j: float):
    """
    Constructs a string representation of a quantum state.

    Args:
        n (int): Principal quantum number.
        l (int): Orbital angular momentum quantum number.
        j (float): Total angular momentum quantum number.

    Returns:
        str: A string representing the quantum state in the format "nlj", where
                'n' is the principal quantum number, 'l' is the orbital symbol, and
                'j' is the total angular momentum quantum number.

    Raises:
        IndexError: If the orbital angular momentum quantum number 'l' is out of range
                    for the defined orbital symbols.
    """
    # Define a dictionary to convert orbital symbols to quantum numbers
    orb_mapping = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'k']
    return str(n) + orb_mapping[l] + str(j)


def is_lin_dip_permitted(state_a:tuple, state_b:tuple):
    """
    Check if a transition between two states is permitted under the leading dipole approximation with linear polarization.
    Args:
        state_a (tuple): Quantum numbers (n, l, j, m) for the initial state.
        state_b (tuple): Quantum numbers (n, l, j, m) for the final state.
    Returns:
        bool: True if the transition is permitted, False otherwise.

    """
    (na, la, ja, ma) = state_a
    (nb, lb, jb, mb) = state_b

    if abs(la-lb) != 1 or abs(ja-jb) > 1 or ma != mb : return False
    
    return True


def compare_transitions(ref_delta: float, trans_df: pd.DataFrame):
    """
    Compare the transition differences between a reference energy delta and all other transitions in a DataFrame.
    Parameters:
    ref_trans (float): The reference energy delta.
    trans_df (pd.DataFrame): A DataFrame containing transition data with a 'delta' column.
    Returns:
    pd.DataFrame: A DataFrame sorted by the absolute differences in 'delta' values between the reference transition and all other transitions.
    """

    ref_delta_abs = abs(ref_delta)
    trans_diff_list = []
    for _, row in trans_df.iterrows():
        trans_diff_list.append(abs(abs(row['delta']) - ref_delta_abs))

    trans_diff_df = pd.DataFrame(trans_diff_list, index=trans_df.index, columns=['trans_diff'])
    return trans_diff_df.sort_values(by='trans_diff')


def wavelength_to_frequency(wavelength_nm : float):
    """
    Convert wavelength in nm to frequency in MHz.
    
    Parameters:
    wavelength_nm (float): Wavelength in nanometers.
    
    Returns:
    float: Frequency in megahertz.
    """

    c_light = 299792458  # Speed of light in m/s
    wavelength_m = wavelength_nm * 1e-9  # Convert nm to meters
    frequency_hz = c_light / wavelength_m  # Frequency in Hz
    frequency_mhz = frequency_hz * 1e-6  # Convert Hz to MHz

    return frequency_mhz


def frequency_to_wavelength(frequency_mhz : float):
    """
    Convert frequency in MHz to wavelength in nm.
    
    Parameters:
    frequency_mhz (float): Frequency in megahertz.
    
    Returns:
    float: Wavelength in nanometers.
    """
    
    c_light = 299792458  # Speed of light in m/s
    frequency_hz = frequency_mhz * 1e6  # Convert MHz to Hz
    wavelength_m = c_light / frequency_hz  # Wavelength in meters
    wavelength_nm = wavelength_m * 1e9  # Convert meters to nanometers
    
    return wavelength_nm