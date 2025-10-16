import pycatima as catima
import pandas as pd

"""
Example script demonstrating the use of PyCatima to calculate energy loss of ions in a target material.
Info taken from documentation @ https://github.com/hrosiak/pycatima/blob/main/docs/pycatima.md
"""


# --- User-defined inputs ---
Z_p = 3.0        # projectile Z
A_p = 6.0        # projectile A
E_p = 5.33       # energy in MeV/u
charge_state = 3.0  # charge state

Z_t = 4.0
A_t = 9.0
target_material = [A_t, Z_t, 1.0]  #
thickness = 9.25e-5         # g/cm^2
# -------------------------

"""
Projectile is defined by pycatima.Projectile class. It is initialized using:

Projectile(A, Z, Q=Z, T=0)

1) A is mass in u units
2) Z is proton number
3) Q is charge state
4) T is energy in Mev/u units.

"""
proj = catima.Projectile(A=A_p, Z=Z_p, Q=charge_state, T=E_p)


"""
Material is defined by pycatima.Material class. The recommended way of initialization is usign the following init signature:

Material(elements, density, thickness, i_potential, mass)

    elements - list of elements, where element is defined as list of [A, Z, STN]. 
        A is atomic mass of the element, if 0 natural abundance atomic mass is taken 
        Z is the proton number of the element 
        STN is the stoichiometric if >=1.0 or weight fraction if < 1.0
    density - optional, defaults to 0
    thickness - optional, if not defined 0
    i_potential - optional, if <=0 it will be calulated using Bragg rule from elemental ionization potentials
    mass - optional, if <=0 it will be calculated from elements masses and STN number.

"""
# beryll = catima.Material(material=target_material, thickness=thickness)
# beryll = catima.Material(material=target_material, density=1.85)
# beryll.thickness(thickness)

beryll = catima.get_material(4)
beryll.thickness(thickness)

e_step = catima.dedx(proj, beryll)
print(f"Initial dE/dx: {e_step:.3f} MeV/cm")



# --- Calculate stopping results ---
results = pd.DataFrame({catima.calculate(proj, beryll)})

# --- Inspect Results object ---
print("=== Results from PyCatima calculate() ===")
# print(f"Initial energy: {results.initial_energy} MeV/u")
# print(f"Residual energy: {results.residual_energy} MeV/u")
# print(f"Total energy loss: {results.energy_loss_total} MeV/u")
# print(f"Range: {results.range} g/cm^2")


# Save to CSV
results.to_csv("energy_loss_curve.csv", index=False)
print("Saved energy loss data to energy_loss_curve.csv")

# Optional: print detailed stopping curve
# print("\n--- dE/dx Curve ---")
# for x, dedx in zip(results.depths, results.dedx_curve):
#     print(f"{x:.3f} g/cm^2 -> {dedx:.3f} MeV/(g/cm^2)")
