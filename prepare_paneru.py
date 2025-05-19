import os
import numpy as np

data = np.loadtxt("data/Paneru.dat")

# If AZURE2 directory does not exist, create it
if not os.path.exists('data/SAMMY'):
    os.makedirs('data/SAMMY')

# Multiply the energy by 1e6
data[:,0] *= 1e6

# Correct projectile and ejectile names
projectile = "h"
ejectile = "h"

# Round the angle to 2 decimal places
for i in range(len(data)):
    data[i,1] = round(data[i,1], 2)

# Now we have to create one file for each angle in the data (second column)
for angle in np.unique(data[:,1]):
    # Get the rows with the same angle
    rows = np.where(data[:,1] == angle)[0]
    # Get the data for the angle
    data_angle = data[rows]
    # Angle in string format
    angle_str = '{:.2f}'.format(angle)
    # Create the name for the file
    name = 'Paneru_' + projectile + ejectile + '_' + angle_str + '.twenty'
    print(name)
    # The data must be saved with a line with energy and cross section
    # Then a second line with tab and the error
    # Then the second data point and so on
    # Open the file
    with open('data/SAMMY/' + name, 'w') as f:
        for row in data_angle:
            f.write('{:d}. {:.3e}\n'.format(int(row[0]), row[2]))
            f.write('\t\t {:.3e}\n'.format(row[3]))