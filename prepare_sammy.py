import os
import numpy as np

from PoPs import database as databaseModule

defaultPops = './ripl2pops.xml'
pops = databaseModule.Database.read( defaultPops )

def getmz(nucl):
    if nucl[0]=='2':
        multiplicity = 2 
        nucl = nucl[1:]
    else: 
        multiplicity = 1
        
    n = pops[nucl.replace('-','')]
    mass = n.getMass('amu')
    if nucl=='n' or nucl=='photon':
        charge = 0
    else:
        charge = n.nucleus.charge[0].value
    return (mass,charge)

def tolab(m1, m2, e):
    return (m1+m2)/m2 * e

def tocm(m1, m2, e):
    return m2/(m1+m2) * e

def rutherford(z1, z2, energy, angle):
    angle = angle * np.pi / 180
    k = 1.296 * pow( z1 * z2 / energy, 2 )
    ruth = k / pow( np.sin( angle / 2 ), 4 )
    return ruth

def cm_to_lab_angle(theta_cm_deg, m1, m2):
    theta_cm_rad = np.radians(theta_cm_deg)
    num = np.sin(theta_cm_rad)
    den = np.cos(theta_cm_rad) + m1 / m2
    theta_lab_rad = np.arctan2(num, den)
    return np.degrees(theta_lab_rad)

def lab_to_cm_angle(theta_lab_deg, m1, m2):
    theta_lab_rad = np.radians(theta_lab_deg)
    num = np.sin(theta_lab_rad)
    den = np.cos(theta_lab_rad) - m1 / m2
    theta_cm_rad = np.arctan2(num, den)
    return np.degrees(theta_cm_rad)

# Define the path to the data
path = 'data/X4'

# Read the datafile.props.csv with numpy, where comma is the delimiter, skip the header
props = np.genfromtxt(path + '/datafile.props.csv', delimiter=',', dtype=str, skip_header=1)

# Get the list of files in the folder
files = os.listdir(path)

# Loop over the files
for file in files:

    # If not .dat file, skip
    if not file.endswith('.dat'):
        continue

    # Check the row in the props file by looking at the fifth column which is the filename
    row = np.where(props[:,5] == file)[0]
    metadata = props[row][0]

    is_cm = True if metadata[17] == 'cm' else False
    is_ruth = True if metadata[21] == 'TRUE' else False
    is_lab_angle = True if metadata[12] == 'TRUE' else False
    is_ang_integrated = True if metadata[8] == 'TRUE' else False

    m1, z1 = getmz(metadata[0])
    m2, z2 = getmz(metadata[2])

    print( 'Preparing file ' + file + ' for AZURE2')

    # Read the data
    data = np.loadtxt(path + '/' + file)

    # If only one row is present, we get a (4,) array instead of (1,4)
    if len(data.shape) == 1:
        data = np.array([data])

    # If center of mass system, convert the energy to lab system
    if is_cm:
        data[:,0] = tocm(m1, m2, data[:,0])

    # If ratio to Rutherford, convert the cross section
    if( is_ruth ):
        for i in range(len(data)):
            energy = tocm(m1, m2, data[i,0])
            angle = data[i,1]
            conversion = rutherford(z1, z2, energy, angle)
            data[i,2] = data[i,2] * conversion
            data[i,3] = data[i,3] * conversion

    # If the angle is in CM, convert it to lab angle
    if not is_lab_angle and not is_ang_integrated:
        angle_original = data[:,1].copy()
        for i in range(len(data)):
            angle = data[i,1].copy()
            cross = data[i,2].copy()
            cross_error = data[i,3].copy()
            # Ivanovic uses 3He instead of 4He as target, so we convert the energies
            if( "A1014010" in file ):
                angle_lab = cm_to_lab_angle(angle, m2, m1)
                # Now we need to convert the cross section to lab system, since the infinitesimal is different
                gamma = m2 / m1 # Iliadis "Nuclear Physics of Stars", Eq. C39
                conversion = ( 1 + gamma * np.cos( np.radians( angle ) ) ) / pow( 1 + gamma**2 + 2 * gamma * np.cos( np.radians( angle ) ), 1.5 ) # Iliadis "Nuclear Physics of Stars", Eq. C44
                cross_lab = cross / conversion
                cross_err_lab = cross_error / conversion
            else:
                angle_lab = cm_to_lab_angle(angle, m1, m2)
                # Now we need to convert the cross section to lab system, since the infinitesimal is different
                gamma = m1 / m2 # Iliadis "Nuclear Physics of Stars", Eq. C39
                conversion = ( 1 + gamma * np.cos( np.radians( angle ) ) ) / pow( 1 + gamma**2 + 2 * gamma * np.cos( np.radians( angle ) ), 1.5 ) # Iliadis "Nuclear Physics of Stars", Eq. C44
                cross_lab = cross / conversion
                cross_err_lab = cross_error / conversion
            data[i,1] = angle_lab
            data[i,2] = cross_lab
            data[i,3] = cross_err_lab

    # Subsitute all -1 in the second column by 0
    data[:,1] = np.where(data[:,1] == -1, 0, data[:,1])

    # If AZURE2 directory does not exist, create it
    if not os.path.exists('data/SAMMY'):
        os.makedirs('data/SAMMY')

    #  Cut the energies below 0.12 MeV for Abramovich (due to screening)
    if( "A0244002" in file ):
        data = data[data[:,0] > 0.119]

    #  Cut the energies above 7 MeV for Harrison (difficult to fit)
    if( "C1003002" in file ):
        data = data[data[:,0] < 7.1]

    # Ivanovic uses 3He instead of 4He as target, so we convert the energies
    if( "A1014010" in file ):
        data[:,0] = data[:,0] * 3.01603 / 4.0026

    # Multiply the energy by 1e6
    data[:,0] *= 1e6

    # Multiply the cross section by 1e-3
    data[:,2] *= 1e-3
    data[:,3] *= 1e-3

    # Where errors are missing, consider 10% error
    data[:,3] = np.where(data[:,3] == 0, data[:,2] * 0.1, data[:,3])

    # Split the publication name
    publication = file.split('-')[0]
    projectile = metadata[0]
    ejectile = metadata[1]
    target = metadata[2]

    # Ivanovic uses 3He instead of 4He as target, so we need to change the projectile
    if( "A1014010" in file ):
        projectile = 'He3'
        ejectile = 'He4'

    # Correct projectile and ejectile names
    if projectile == "photon": projectile = "G"
    if ejectile == "photon": ejectile = "G"
    if projectile == "He3": projectile = "h"
    if ejectile == "He3": ejectile = "h"
    if projectile == "He4": projectile = "a"
    if ejectile == "He4": ejectile = "a"
    if projectile == "H1": projectile = "p"
    if ejectile == "H1": ejectile = "p"

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
        name = publication + '_' + projectile + ejectile + '_' + angle_str + '.twenty'

        print(name)

        # The data must be saved with a line with energy and cross section
        # Then a second line with tab and the error
        # Then the second data point and so on

        # Open the file
        with open('data/SAMMY/' + name, 'w') as f:
            for row in data_angle:
                f.write('{:d}. {:.3e}\n'.format(int(row[0]), row[2]))
                f.write('\t\t {:.3e}\n'.format(row[3]))