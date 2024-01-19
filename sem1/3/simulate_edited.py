"""
Symplectic Euler and Velocity Verlet time integration of
a particle moving in a double well potential.

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

The potential is V(x) = a*x^4 - b*x^2, where
a and b are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.

Author: Paras Ladwa
Number: s2188899

"""

import sys
import math
import scipy
import numpy as np
import matplotlib.pyplot as pyplot
#from particle1D import Particle1D
#import 3d
from particle3D import Particle3D



def pull(i, particle_file):
    '''
    Parameters
    ----------
    i : integer
        value of which line of the file to read
    particle_file : str
        which file to open and read from

    Returns ith line from particle_file
    -------

    '''  
    file = open(particle_file, 'r')
    f = file.readlines()
    return f[i]
    
    file.close()




def make_particle(line, particle_file):
    '''

    Parameters
    ----------
    line : int
        line of file to use
    particle_file : str
        file name of particle data.

    Returns instance of particle with data from file, using particle3D
    '''
    data = pull(line, particle_file)
    return Particle3D.read_line(data)



def extract_vectors(p1, p2):
    '''

    Parameters
    ----------
    p1 : Particle3D
    p2 : Particle3D

    Returns
    -------
    r_1_vector : position vector
    r_2_vector : ''
    r_12_vector : relative posn
    r_12_mag : magnitude of relative posn
    r_12_hat : unit vect of rel posn
    '''
    
    #extract position vectors
    r_1_vector , r_2_vector = p1.position, p2.position
    
    #relative position vector
    r_12_vector = r_2_vector - r_1_vector
    
    #magnitude of r_12
    r_12_mag = np.linalg.norm(r_12_vector)
    
    #unit vector of r_12_vector
    r_12_hat = r_12_vector / r_12_mag
    
    return r_1_vector, r_2_vector, r_12_vector, r_12_mag, r_12_hat



def force_double_well(p1, p2, alpha, D_e, r_e):
    """
    Return the force on a particle in a Morse potential.

    The force is given by
        F1(r1, r2) = 2 alpha D_e [ 1 - exp-alpha(r_12 - r_e) ] exp-alpha(r_12 - r_e)*r_hat

    Parameters
    ----------
    p1: Particle3D
    p2: Particle3D
    alpha: float
    D_e: float
    r_e: float

    Returns
    -------
    force: vector
    """
    
    r_12_mag, r_12_hat = extract_vectors(p1, p2)[-2:]
    
    #computing force using given eqn
    force = 2*alpha*D_e*(1-math.exp(-1*alpha*(r_12_mag-r_e)))*(math.exp(-alpha*(r_12_mag-r_e))*r_12_hat)
    

    return force


def potential_double_well(p1, p2, alpha, D_e, r_e):
    """
    Method to return morse potential energy
    of particles

    Parameters
    -----------
    p1: Particle3D
    p2: Particle3D
    alpha: float
    D_e: float
    r_e: float

    Returns
    -------
    potential: float
    """
    
    r_12_mag = extract_vectors(p1, p2)[3]
    
    potential = D_e*(((1-math.exp(-1*alpha*(r_12_mag-r_e)))**2)-1)
    return potential



def average_difference_peaks(input_list):
    '''

    Parameters
    ----------
    input_list : list
        list of data

    Returns
    -------
    wavenumber : float
        wave number

    '''
    #uses scipy to find maxima
    indicies = scipy.signal.find_peaks(input_list)
    
    differences = []
    
    #iterate through list of maxima, finding the wavelength
    #of each oscillation to find the average
    for i in range(0, len(indicies[0])-1):
        differences.append(indicies[0][i+1] - indicies[0][i])
    
    #converting to wavenumber in correct units
    average = (sum(differences)/len(differences))*(1.018050571*(10**-14))
    frequency = 1/average
    wavenumber = frequency/(3*(10**8))
    
    return wavenumber


def data_inaccuracy(initial, data):
    '''

    Parameters
    ----------
    initial : float
        initial energy.
    data : list 
        lsit of energies during oscillations
        

    Returns
    -------
    inaccuracy : float
        E_0 / delta E.

    '''
    data_range = max(data) - min(data)
    inaccuracy = data_range/initial
    return inaccuracy
    

'''
same code but replaces sysargvs as parameters
and returns necessary vals
so i can use it in loops within test_edited.py

'''
def main(mode, particle_file, outfile_name, dt):
    '''
    Parameters
    ----------
    mode : str
        which integration scheme to use
    particle_file : str
        particle file name
    outfile_name : str
        output file name.
    dt : float
        value of dt to step through each iteration.

    Returns
    -------
    wavenumber - float
        wavenumber.
    energy uncertainty - float
        energy uncertainty
    energies : 1d nparray
        list of energies

    '''


    dt = np.float64(dt)
        
    # Open the output file for writing ("w")
    outfile = open(outfile_name, "w")

    # Choose our simulation parameters
    
    consts = pull(0, particle_file)
    consts = consts.split()
    D_e, r_e, alpha = float(consts[0]), float(consts[1]), float(consts[2])
    
    
    #using a dictionary for number of iterations needed
    #for each loop of simulation
    numstep_dict = {
        0.5 : 10000,
        0.1 : 100000,
        0.01 : 100000,
        0.001 : 100000,
        0.0001 : 100000,
        0.00001 : 1100000
        }
    numstep = numstep_dict[dt]
    
    #feed_find_peaks - to feed every nth datapoint to
    #average_difference_peaks() function
    #otherwise find_peaks seems to not work
    #unsure why
    feed_find_peaks = int(0.01 / dt)
    
    #for cases where feed_find_peaks == 0
    if feed_find_peaks == 0:
        feed_find_peaks = 1
    
    #initialise time val
    time = 0
    

    
    # Set up particle initial conditions:
    #  position x0 = 0.0
    #  velocity v0 = 1.0
    #  mass      m = 1.0
    
    p1 = make_particle(1, particle_file)
    p2 = make_particle(2, particle_file)


    # Get initial force
    force = force_double_well(p1, p2, alpha, D_e, r_e)

    # Write out starting time, position, and energy values
    # to the output file.
    energy = p1.kinetic_energy() + p2.kinetic_energy() + potential_double_well(p1, p2, alpha, D_e, r_e)
    initial_energy = energy
    outfile.write(f"{time}    {extract_vectors(p1, p2)[3]}    {energy}\n")

    # Initialise numpy arrays that we will plot later, to record
    # the trajectories of the particles.
    times = np.zeros(numstep)
    
    #positions = np.zeros(numstep)
    positions = np.zeros(numstep) #changed for correct dimensionality when adding at end of loop
    energies = np.zeros(numstep)

    # Start the time integration loop
    for i in range(numstep):
        # Update the positions and velocities.
        # This will depend on whether we are doing an Euler or verlet integration
        if mode == "euler":
            '''
            my added notes
             - update posn using vel 
             - using new posns find force thus vels
            '''
            # Update particle position
            p1.update_position_1st(dt)
            p2.update_position_1st(dt)
            
            # Calculate force
            force = force_double_well(p1, p2, alpha, D_e, r_e)
            

            # Update particle velocity 
            p1.update_velocity(dt, force)
            p2.update_velocity(dt, -force)

        elif mode == "verlet":
            '''
            veloctiy update is different
            '''
            # Update particle position using previous force
            p1.update_position_2nd(dt, force)
            p2.update_position_2nd(dt, -force)
            
            # Get the force value for the new positions
            force_new = force_double_well(p1, p2, alpha, D_e, r_e)

            # Update particle velocity by averaging
            # current and new forces
            
            p1.update_velocity(dt, ((0.5)*(force+force_new)))
            p2.update_velocity(dt, (-(0.5)*(force+force_new)))
            
            # Re-define force value for the next iteration
            force = force_new
        else:
            raise ValueError(f"Unknown mode {mode} - should be euler or verlet")

        # Increase time
        time += dt
        
        # Output particle information
        energy = p1.kinetic_energy() + p2.kinetic_energy() + potential_double_well(p1, p2, alpha, D_e, r_e)
        outfile.write(f"{time} {extract_vectors(p1, p2)[3]} {energy}\n")

        # Store the things we want to plot later in our arrays
        times[i] = time
        positions[i] = extract_vectors(p1, p2)[3]
        energies[i] = energy

    # Now the simulation has finished we can close our output file
    outfile.close()

    # Plot particle trajectory to screen. There are no units
    # here because it is an abstract simulation, but you should
    # include units in your plot labels!
    pyplot.figure()
    pyplot.title(f"{mode} {particle_file[:-4]}: seperation vs time")
    pyplot.xlabel('Time')
    pyplot.ylabel('seperation')
    pyplot.scatter(times[::feed_find_peaks], positions[::feed_find_peaks])
    pyplot.show()
    pyplot.close()

    # Plot particle energy to screen
    pyplot.figure()
    pyplot.title(f"{mode} {particle_file[:-4]}: total energy vs time")
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(times, energies)
    pyplot.show()


    #return necessary values
    return data_inaccuracy(initial_energy, energies), average_difference_peaks(positions[::feed_find_peaks]), energies





















