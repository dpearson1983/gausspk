import numpy as np

# Define a data type to store complex numbers for ease
dt = np.dtype([('re', np.float64), ('im', np.float64)])

# Load the binary file
dk = np.fromfile('dk_0001.bin', dtype=dt)

# Print one of the complex numbers to demonstrate how to access
print("(" + str(dk[1]['re']) + ", " + str(dk[1]['im']) + ")\n")
