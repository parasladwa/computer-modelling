"""
TODO: approximate partial solar system simulation
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
import sys

#import basic functions
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




def particles_from_file(filename = "mini_system.txt"):
    '''
    Parameters
    ----------
    filename : str, opt
        filename of particles.
        The default is "mini_system.txt".

    Returns
    -------
    particles : class
        list of instances of particles, Particle3D

    '''
    
    particles = []
    
    file = open(filename, 'r')
    f = file.readlines()
    
    for particle in f:
        particles.append(Particle3D.read_line(particle))
    
    return particles
        
    


    
    
def main():
    '''
    implementing verlet velocity integration scheme
    to reduced solar system, running simulation of a year,
    plotting appropriate positions and energies.

    '''
    # Set up simulation parameters:
    # about 10 years of simulation with 1 day timesteps
    
    # numstep dt inputfile outputxyz 
    # %run template.py 36500 0.01 mini_system.txt outfile.xyz
    
    if len(sys.argv) in [5, 6]:
        
        numstep = int(sys.argv[1])
        dt = float(sys.argv[2])
        particle_file = sys.argv[3]
        outfile_name = sys.argv[4]
        
        
        if len(sys.argv) == 6:
            extra_out = sys.argv(5)
    
        print("SYS ARGVS OK")
    else:
        print("Incorrect sys argv's\nCorrect form as follows :")
        print("%run template.py <numstep> <dt> <particle file> <xyz outfile> <OPTIONAL out>")
        sys.exit(1)
        



    
    #open outfile
    outfile = open(outfile_name, 'w')
    
    # Initial conditions of the system
    particles = particles_from_file(particle_file)
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

    #compute initial forces for first loop iteration
    separations = b_f.compute_separations(particles)
    forces, potential = b_f.compute_forces_potential(particles, separations)

    # Main time integration loop
    for i in range(numstep):
        times[i] = time
        time += dt
        
        #outfile headers
        if i%250 == 0:
            outfile.write(f"{n}\n")
            outfile.write(f"Point = {i}\n")
            
            
            for particle in particles:
                outfile.write(f"{particle.label} {str(particle.position)[1:-1]}\n")
        
        
        #update all particle positions
        #store particle positions in array
        for j, particle in enumerate(particles):
            particle.update_position_2nd(dt, forces[j])
            positions[j][i] = particle.position #maybe error here
            
         
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


    
    '''
    SORT OUT PLOTS
    '''
    
    
    
    #dictionary of particle locations
    particle_dict = {}
    for i, particle in enumerate(particles):
        particle_dict[particle.label] = i
        particle_dict[i] = particle.label
    print(particle_dict)
    
    
    
    
    
    def get_positions(particle_label):
        '''
        Parameters
        ----------
        particle_label : str
            name of particle (particle.label)

        Returns
        -------
        data : np array
            list of positions of particle throughout time
            such that data[x] is all the x values of the 
            position throughout the simulation
        
        '''
        index = particle_dict[particle_label]
        data = positions[index, :]  
        data = np.transpose(data)
        x, y, z = data
        return x, y, z

  

    
    # Make two plots, of the Mercury - Sun x distance,
    # and of the trajectory x-vs-y of the Earth's orbit.
    pyplot.title('Mercury-Sun Location')
    pyplot.xlabel('time / days')
    pyplot.ylabel('x / AU')
    pyplot.plot(times, get_positions("Mercury")[0] - get_positions("Sun")[0])
    pyplot.show()


    pyplot.title('Earth Trajectory')
    pyplot.xlabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(get_positions("Earth")[0],  get_positions("Earth")[1])
    pyplot.show()


    pyplot.title('Total Energy')
    pyplot.xlabel('x / days')
    pyplot.ylabel('Energy / M_earth AU^2 / day^2')
    pyplot.plot(times, energy)
    pyplot.show()


    #add a plot of the position of the sun with time
    pyplot.title('Sun Trajectory')
    pyplot.ylabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(get_positions("Sun")[0], get_positions("Sun")[1])
    pyplot.show()


    #add a plot of the orbit of the trajectory of the moon around the earth
    pyplot.title('Moon - Earth Location')
    pyplot.xlabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(get_positions("Moon")[0] - get_positions("Earth")[0], get_positions("Moon")[1] - get_positions("Earth")[1])
    pyplot.show()
    

    # You can add other useful plots here to check the system.
    
    
    
    
    
    
    
    
    
    
    #UNIT 4 WORK
    def energy_deviation(energies):
        initial = energies[0]
        minimum = min(energies)
        maximum = max(energies)
        
        deviation = abs((maximum - minimum)/initial)
        
        return deviation
    
    print(f"Energy deviation = {energy_deviation(energy)}")
    
    
    
    
    
    def helions(big, small):
        perihelion = np.inf #minimum
        aphelion = 0 #maximum
        
        pos_big = get_positions(big)
        pos_small = get_positions(small)
        
        perihelion = 0
        aphelion = np.inf
        
        for i in range(0, len(pos_big)):
            
            distance = abs(pos_big[i] - pos_small[i])
            
            if distance < perihelion:
                perihelion = distance
                print(i)
            
            if distance > aphelion:
                aphelion = distance
                print(i)
        return perihelion, aphelion
    
    print(helions("Sun", "Earth"))
    
    # def helion_lists():
        
    #     orbit_sun = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
    
        
    #     for p in particles:
    #         if p not in orbit_sun:
    #             continue
    #         perihelion = 
    #         #here
        
        
    
    
    outfile.close()
    

    



# This python standard code makes it so that the "main"
# function is only called when you run the file directly,
# not when you just import it from another python file.
if __name__ == "__main__":
    main()

