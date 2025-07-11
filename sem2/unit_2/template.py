"""
approximate partial solar system simulation
using verlet integration scheme,
partial due to smaller number of planets, listed in function below
plots certain orbits using matplotlib.

uses particle3D class where each instance is a planet

Paras Ladwa
s2188899
"""

import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D

# importing basic functions
import basic_functions as b_f


def generate_simple_solar_system():
    # This is a smaller test set of particles to use in tests.
    particles = [
        Particle3D('Sun', 332946.0, np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])),
        Particle3D('Mercury', 0.055, np.array([0.08051955462645541, -0.4480662860374615, -0.04400182149190969 ]), np.array([0.02204825304859444, 0.00641708881027627, -0.00149796026125163])),
        Particle3D('Earth', 1.0, np.array([-0.4858557302072681, -0.8881256929389151, 4.644518255972374e-05]), np.array([0.01482059753079534, -0.008322015184828249, 8.879147180855189e-08])),
        Particle3D('Moon', 0.0123, np.array([-0.4863273339804307, -0.8855079737558565, 0.0002690858287502895]), np.array([0.01425788196666878, -0.008404177858858091, 2.210145855707591e-05]) ),
    ]
    return particles

    
    
def main():
    '''
    implementing verlet velocity integration scheme
    to reduced solar system, running simulation of a year,
    plotting appropriate positions and energies.

    '''

    
    # Set up simulation parameters:
    # about 10 years of simulation with 1 day timesteps
    dt = 0.01  
    numstep = 36500

    # Initial conditions of the system
    particles = generate_simple_solar_system()
    time = 0.0

    #subtract the centre-of-mass velocity
    com_vel = Particle3D.com_velocity(particles)
    for particle in particles:
        particle.velocity -= com_vel

    # Initialise arrays that we will store results in
    n = len(particles)
    times = np.zeros(numstep)
    energy = np.zeros(numstep)
    positions = np.zeros((n, numstep, 3))

    # TODO done: compute initial forces for first loop iteration
    separations = b_f.compute_separations(particles)
    forces, potential = b_f.compute_forces_potential(particles, separations)

    # Main time integration loop
    for i in range(numstep):
        times[i] = time
        time += dt
        
        #update all particle positions
        #store particle positions in array
        for j, particle in enumerate(particles):
            particle.update_position_2nd(dt, forces[j])
            positions[j][i] = particle.position
            
        
        #get new separations and new forces on all particles, and the potential
        separations = b_f.compute_separations(particles)
        forces, potential = b_f.compute_forces_potential(particles, separations)
        
        
        #update all particle velocities
        for k, particle in enumerate(particles):
            particle.update_velocity(dt, forces[k])
            
        
        
        #replace forces with new forces for next iteration
        separations = b_f.compute_separations(particles)
        forces, potential = b_f.compute_forces_potential(particles, separations)
    
        #compute the kinetic energy and save the total energy
        energy[i] = Particle3D.total_kinetic_energy(particles)

        
    # Make two plots, of the Mercury - Sun x distance,
    # and of the trajectory x-vs-y of the Earth's orbit.
    pyplot.title('Mercury-Sun Location')
    pyplot.xlabel('time / days')
    pyplot.ylabel('x / AU')
    pyplot.plot(times, positions[1, :, 0] - positions[0, :, 0])
    pyplot.show()

    pyplot.title('Earth Trajectory')
    pyplot.xlabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(positions[2, :, 0],  positions[2, :, 1])
    pyplot.show()

    pyplot.title('Total Energy')
    pyplot.xlabel('x / days')
    pyplot.ylabel('Energy / M_earth AU^2 / day^2')
    pyplot.plot(times, energy)
    pyplot.show()


    #plot of the position of the sun with time
    pyplot.title('Sun Trajectory')
    pyplot.ylabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(positions[0, :, 0], positions[0, :, 1])
    pyplot.show()

    #plot of the orbit of the trajectory of the moon around the earth
    pyplot.title('Moon - Earth Location')
    pyplot.xlabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(positions[3, :, 0] - positions[2, :, 0], positions[3, :, 1] - positions[2, :, 1])
    pyplot.show()

    

# This python standard code makes it so that the "main"
# function is only called when you run the file directly,
# not when you just import it from another python file.
if __name__ == "__main__":
    main()

