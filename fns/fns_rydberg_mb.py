import numpy as np
import scipy as sp
import qutip as qtp
import re

def parse_spec_state(state_str):
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


def build_state_str(n:int, l:int, j:float):
    # Define a dictionary to convert orbital symbols to quantum numbers
    orb_mapping = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'k']
    return str(n) + orb_mapping[l] + str(j)


def is_lin_dip_permitted(state_a:list, state_b:list):
    """
    Boolean check if a transition is permitted (not forbidden)
     at leading dipole approx under linear polarization.

    """
    [na, la, ja, ma] = state_a
    [nb, lb, jb, mb] = state_b

    if abs(la-lb) != 1 or abs(ja-jb) > 1 or ma != mb : return False
    
    return True
