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
    return (m1+m2)/m1 * e

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

def transform_dsigma_cm_to_lab(theta_cm_deg, dsigma_domega_cm, m1, m2):
    theta_lab_deg = cm_to_lab_angle(theta_cm_deg, m1, m2)
    theta_lab_rad = np.radians(theta_lab_deg)
    theta_cm_rad = np.radians(theta_cm_deg)

    # Derivative d(theta_lab)/d(theta_cm) using small angle difference
    delta = 1e-5  # small angle in degrees
    dtheta_cm1 = theta_cm_deg - delta
    dtheta_cm2 = theta_cm_deg + delta

    theta_lab1 = cm_to_lab_angle(dtheta_cm1, m1, m2)
    theta_lab2 = cm_to_lab_angle(dtheta_cm2, m1, m2)

    dtheta_lab_dtheta_cm = (theta_lab2 - theta_lab1) / (2 * delta)

    jacobian = np.abs(
        (np.sin(theta_cm_rad) / np.sin(theta_lab_rad)) * dtheta_lab_dtheta_cm
    )

    dsigma_domega_lab = dsigma_domega_cm * jacobian

    return dsigma_domega_lab

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
        data[:,0] = tolab(m2, m1, data[:,0])

    # If ratio to Rutherford, convert the cross section
    if( is_ruth ):
        for i in range(len(data)):
            energy = tocm(m1, m2, data[i,0])
            angle = data[i,1]
            conversion = rutherford(z1, z2, energy, angle)
            data[i,2] = data[i,2] * conversion
            data[i,3] = data[i,3] * conversion

    # Subsitute all -1 in the second column by 0
    data[:,1] = np.where(data[:,1] == -1, 0, data[:,1])

    # If AZURE2 directory does not exist, create it
    if not os.path.exists('data/AZURE2'):
        os.makedirs('data/AZURE2')

    # Multiply the cross section by 1e-3
    data[:,2] *= 1e-3
    data[:,3] *= 1e-3

    # Where errors are missing, consider 10% error
    data[:,3] = np.where(data[:,3] == 0, data[:,2] * 0.1, data[:,3])

    # Now sort the data by the first column first
    data = data[data[:,1].argsort()]

    #  Cut the energies below 0.12 MeV for Abramovich (due to screening)
    if( "A0244002" in file ):
        data = data[data[:,0] > 0.119]

    #  Cut the energies above 7 MeV for Harrison (difficult to fit)
    if( "C1003002" in file ):
        data = data[data[:,0] < 7.1]

    # Ivanovic uses 3He instead of 4He as target, so we convert the energies
    if( "A1014010" in file ):
        data[:,0] = data[:,0] * 3.01603 / 4.0026

    # Save the data with constatnt exponential notation with fixed precision
    np.savetxt('data/AZURE2/' + file, data, fmt='%.3e', delimiter=' ')