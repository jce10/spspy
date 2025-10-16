from NuclearData import get_excitations
from NuclearData import generate_nucleus_id

print("testing 13C ")
print(generate_nucleus_id(6, 13))



# Test cases
print("Testing 13C...")
print(get_excitations(6, 13))  # Carbon-13 (Z=6, A=13)

print("Testing 12C...")
print(get_excitations(6, 12))  # Carbon-12

print("Testing 6Li...")
print(get_excitations(3, 6))   # Lithium-6