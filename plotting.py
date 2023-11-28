import matplotlib.pyplot as plt
import numpy as np

# Read data from the text file
data = []
with open('data.txt', 'r') as file:
    for line in file:
        x, y, z = map(float, line.split())
        data.append((x, y, z))

# Separate data into three arrays
data = np.array(data)
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

plt.figure(figsize=(8, 6))
plt.plot(x, y, color='b', label="Mode 1")
plt.plot(x, z, color='r', label="Mode 2")
plt.xlabel('Mass of Each Particle (kg)', fontsize=14)
plt.ylabel('Angular Frequency (s^-1)', fontsize=14)


plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.legend(fontsize=12)
plt.grid(True)

plt.show()
