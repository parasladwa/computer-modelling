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
import math
import scipy
#import scipy.optimize
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
    
        print("\n   SYS ARGVS OK")
        
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
        if i%2500 == 0:
            outfile.write(f"{n}\n")
            outfile.write(f"Point = {i}\n")
            
            
            for particle in particles:
                outfile.write(f"{particle.label} {str(particle.position)[1:-1]}\n")
        
        #update all particle positions
        #store particle positions in array
        for j, particle in enumerate(particles):
            particle.update_position_2nd(dt, forces[j])
            positions[j][i] = particle.position #maybe error here
        
        previous_force = forces
            
        #get new separations and new forces on all particles, and the potential
        separations = b_f.compute_separations(particles)
        forces, potential = b_f.compute_forces_potential(particles, separations)
        
        
        #update all particle velocities
        average_force = (previous_force + forces)/2
        for k, particle in enumerate(particles):
            particle.update_velocity(dt, average_force[k])
            
        #replace forces with new forces for next iteration
        separations = b_f.compute_separations(particles)
        forces, potential = b_f.compute_forces_potential(particles, separations)
    
        #compute the kinetic energy and save the total energy
        energy[i] = Particle3D.total_kinetic_energy(particles) + potential


    
    '''
    SORT OUT PLOTS
    '''
        
    #dictionary of particle locations
    #maybe move this into get_positions function to 
    #keep local, unless used elsewhere later
    particle_dict = {}
    for i, particle in enumerate(particles):
        particle_dict[particle.label] = i
        particle_dict[i] = particle.label
    
    

    
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
    # pyplot.title('Mercury-Sun Location')
    # pyplot.xlabel('time / days')
    # pyplot.ylabel('x / AU')
    # pyplot.plot(times, get_positions("Mercury")[0] - get_positions("Sun")[0])
    # pyplot.show()


    # pyplot.title('Earth Trajectory')
    # pyplot.xlabel('x / AU')
    # pyplot.ylabel('y / AU')
    # pyplot.plot(get_positions("Earth")[0],  get_positions("Earth")[1])
    # pyplot.show()


    # pyplot.title('Total Energy')
    # pyplot.xlabel('x / days')
    # pyplot.ylabel('Energy / M_earth AU^2 / day^2')
    # pyplot.plot(times, energy)
    # pyplot.show()


    # #add a plot of the position of the sun with time
    # pyplot.title('Sun Trajectory')
    # pyplot.ylabel('x / AU')
    # pyplot.ylabel('y / AU')
    # pyplot.plot(get_positions("Sun")[0], get_positions("Sun")[1])
    # pyplot.show()


    # #add a plot of the orbit of the trajectory of the moon around the earth
    # pyplot.title('Moon - Earth Location')
    # pyplot.xlabel('x / AU')
    # pyplot.ylabel('y / AU')
    # pyplot.plot(get_positions("Moon")[0] - get_positions("Earth")[0], get_positions("Moon")[1] - get_positions("Earth")[1])
    # pyplot.show()
    

    # You can add other useful plots here to check the system.
    
    
    
    
    
    
    
    
    
    
    #UNIT 4 WORK
    def energy_deviation(energies):
        '''
        Parameters
        ----------
        energies : np array
            energy of the system throughout
            the simulation.

        Returns
        -------
        deviation : float
            energy deviation

        '''
        initial = energies[0]
        minimum = min(energies)
        maximum = max(energies)
        
        deviation = abs((maximum - minimum)/initial)
        
        return deviation
    
    print(f"Energy deviation = {energy_deviation(energy)}")
    
    

    
    
            
    def find_central_body(particles, positions):
        """ 
        finds the object closest to the center of the orbits.
        doing this by finding the closest object at a number
        of timesteps and confirming the object is consistent.

        Args:
            particles (list)): list of particle3d instances 
                                (planets, stars etc...)
                                
            positions (np array): full list of particles and 
                                their positions throughout the
                                simulation -------------------------------- dont need this?
                                
        Returns:
            closest[0] (particle3D): object closest to 
                                center of orbits
            
            None (NoneType): if no object is determined to be at the 
                                center of the orbits.
        """
        def closest_body(particles, positions, numstep = 1000):
            positions_at_numstep = []
            
            #list of particles at numstep
            for p in particles:
                current = np.transpose(get_positions(p.label))
                positions_at_numstep.append(current[numstep])
            
            #find center of mass
            mass_position = 0
            net_mass = 0
            for i, pos in enumerate(positions_at_numstep):
                #ith particle corresponds to positions_at_numstep[i]
                mass_position += particles[i].mass*pos
                net_mass += particles[i].mass
            com = mass_position/net_mass
            #print(f"center of mass at {numstep/100} days = {com}")
                
            #find closest particle to com
            closest = [None, np.inf]
            for i, pos in enumerate(positions_at_numstep):
                distance = np.linalg.norm(pos - com)
                if distance < closest[1]:
                    closest = [particles[i], distance]
            #print(f"the closest particle at {numstep/100} days = {closest[0].label}")
            
            return closest[0]
    
        #times to check
        #checking further times as the system
        #will converge to an equilibrium
        times_to_check = [round(num/2) for num in np.linspace(round(numstep/2), numstep-1, 10)]
        closest = []
        for i in times_to_check:
            closest.append(closest_body(particles, positions, i))
        
        
        if len(set(closest)) == 1:
            return closest[0]
        
        else:
            return None
        #come up with messages and what to do with other cases



    
    def apsides():
        """takes all particles and finds min and max 
            distances of orbits from central body.
            also returns pairs of particles which orbit eachother

        Returns:
            particle_pairs (2d list) : list of pairs
                            of particles which orbit
                            eachother
        """

        central_body = find_central_body(particles, positions)
        #failcase
        if central_body == None:
            print(f"No central body was observed\
                    thus no apsides will be computed")
            return 0
        
        particle_pairs = []        
        print(f"\nCentral body identified as {central_body.label}")
        particle_labels = [p.label for p in particles]
        
        def calculate_apsides(p1, p2):
            """ given 2 particles calculates perhelion
                and aphelion manually iterating through
                each timestep
            Args:
                p1 (Particle3D): particle3d instance
                p2 (Particle3D): particle3d instance

            Returns:
                perihelion (float) : minimum distance to central body
                aphelion (float) :   maximum distance to central body
            """
            
            perihelion = np.inf
            aphelion = 0
            posn_1 = np.transpose(get_positions(p1))
            posn_2 = np.transpose(get_positions(p2))
            
            for i in range(0, len(posn_1)):
                
                #calculate distance between pair for each timestep
                distance = np.linalg.norm(posn_1[i] - posn_2[i])
                if distance < perihelion:
                    perihelion = distance
                if distance > aphelion:
                    aphelion = distance
                    
            return perihelion, aphelion
        
        
        #perigee / apogee case
        if  ("Moon" and "Earth") in particle_labels:
            perigee, apogee = calculate_apsides("Moon", "Earth")
            print(f"\nBetween the Moon and Earth :")
            print(f"    Perigee = {perigee} /AU")
            print(f"    Apogee = {apogee} /AU")
            particle_pairs.append(["Moon", "Earth"])
        
        #iterares through above particle pairs
        # calculation orbit maxima and minima
        # using above function and printing to terminal
        for p in particles:
            if (p == central_body) or (p.label == "Moon"):
                continue
            perihelion, aphelion = calculate_apsides(central_body.label, p.label)
            print(f"\nBetween {central_body.label} and {p.label} :")
            print(f"    Perihelion = {perihelion} /AU")
            print(f"    Aphelion = {aphelion} /AU")
            particle_pairs.append([central_body.label, p.label])
            
        return particle_pairs
    
    pairs = apsides()
        
            
    
    
    
    def find_periods():
        
        def orbit_dot_product(pair):
            """
            given a pair of particles calculates dot product
            between initial relative vector and current relative
            vectors at each time. (relative position vectors)
            Args:
                pair (list): 2 instances of particle3d

            Returns:
                x (list) _----------------------------------------------------------------is this needed
                y (1d list of floats) : dot products at each
                                        time as described.
            """
            
            positions_1 = np.transpose(get_positions(pair[0]))
            positions_2 = np.transpose(get_positions(pair[1]))
            initial = positions_1[0] - positions_2[0]
            x, y = [], []
            
            #iterate through positions and log
            #relative dot product
            for i in range(0, len(positions_1)):
                relative = positions_1[i] - positions_2[i]
                y.append(np.dot(initial, relative))
                
            return x, y
        
        
        def period_from_peaks(peaks):
            """given indicies of where peaks in y from above 
                function, calculates period /days

            Args:
                peaks (list of ints): indicies of where maxima
                                        in x is stored

            Returns:
                (float) : averaged
            """
            #failcase
            if len(peaks) < 2:
                return None
            
            #finds average distance between peaks
            #then returns value converted to days
            differences = []
            for i in range(0, len(peaks)-1):
                differences.append(peaks[i+1] - peaks[i])
            average = sum(differences)/len(differences)
            return average*dt
            
        
        #iterates through particle pairs and 
        #prints orbits where applicable to terminal
        for pair in pairs:
            y = orbit_dot_product(pair)[1]
            peaks_indicies = scipy.signal.find_peaks(y)
            period = period_from_peaks(peaks_indicies[0])
            if period == None:
                print(f"\ninsufficient data to deduce orbital perid between {pair[0]} and {pair[1]}")
                continue
            print(f"\norbital period between {pair[0]} and {pair[1]} = {period} days")
            
    find_periods()
    
    
    
    ########
    
    
    # def helions(big, small):
    #     perihelion = np.inf #minimum
    #     aphelion = 0 #maximum
        
    #     pos_big = np.transpose(get_positions(big))
    #     pos_small = np.transpose(get_positions(small))
        
        
    #     for i in range(0, len(pos_big)):
    #         distance = np.linalg.norm(pos_big[i] - pos_small[i])                       
            
    #         if distance < perihelion:
    #             perihelion = distance
            
    #         if distance > aphelion:
    #             aphelion = distance
    #     return perihelion, aphelion
    
    # # p, a = helions("Sun", "Earth")
    # # print(f"{p} -- {a} -- SUN EARTH")
    
    
    # def helion_lists():
        
    #     print("n.b. data may be inaccurate due to incomplete orbits")
    #     def printer(p_1, p_2):
            
    #         min_label, max_label = "Perihelion", "Aphelion"
            
    #         if "Moon" in [p_1, p_2] and "Earth" in [p_1, p_2]:
    #             min_label, max_label = "Perigee", "Apogee"
    #         #min, max
    #         perihelion, aphelion = helions(p_1, p_2)
            
    #         print(f"\n\nBetween {p_1} and {p_2} :")
    #         print(f"    {min_label}   = {perihelion} / AU")
    #         print(f"    {max_label}   = {aphelion} / AU")
   
        
    #     for p in particles:
    #         if p.label != "Sun":
    #             printer("Sun", p.label)
                
        
    #     printer("Earth", "Moon")
        
    #helion_lists()
            
        
    
    
    # def curvefit():
    #     # fail
    #     def orbit_dot_product(p_1, p_2):
            
    #         positions_1 = np.transpose(get_positions(p_1))
    #         positions_2 = np.transpose(get_positions(p_2))
            
    #         initial = positions_1[0] - positions_2[0]
            
    #         x, y = [], []
            
    #         for i in range(0, len(positions_1)):
                
    #             relative = positions_1[i] - positions_2[i]
                
    #             x.append(i/100)
    #             #y.append(math.atan(relative[1] / relative[0]))
                
    #             y.append(np.dot(initial, relative))

    #         return x, y
    #     #CHECK FOR MOON / EARTH
        
    #     """CURVEFIT OPTIMIZATION"""
    #     def sinusiod(x_, amplitude, omega, phi, const):
    #         return amplitude * np.cos(omega * x_ + phi) + const
        
    #     def curve_optimization(x, y):
    #                                                             #[2, 500, 0, 0]
    #         parameters = scipy.optimize.curve_fit(sinusiod, x, y, [2, 0.02, 0, 0])
    #         omega = (parameters[0][1])
    #         T = 2*math.pi /omega
    #         print(f"\nT = {T}")
        
    #         global ynew
    #         ynew = []
    #         for i in x:
    #             ynew.append(sinusiod(i, parameters[0][0], parameters[0][1],\
    #                                     parameters[0][2], parameters[0][3]))        
        
    #     x, y = orbit_dot_product("Earth", "Sun")
    #     curve_optimization(x, y)
        

    #     pyplot.title('dotproduct')
    #     pyplot.xlabel('time / days')
    #     pyplot.ylabel('dot product')
    #     pyplot.plot(x, y, color='blue')
    #     pyplot.plot(x, ynew, color = 'red') 
    #     pyplot.show()
        

    

    



# This python standard code makes it so that the "main"
# function is only called when you run the file directly,
# not when you just import it from another python file.
if __name__ == "__main__":
    main()


# energy units
# label the code
# do i need to pass those params through closest_body ?????????????
# 
# 
# 
# 
# 
# TEST NONE CENTRAL BODY
# 
# 
# comment the code
# 
# 
# why some periods are bad accuracy
# 
# 
# 
# 
# 










################################### UNIT 2 FEEDBACK #######################################
# You haven’t computed a total energy just KE – the magnitude of which does not look right. 

# 2. Correct centre-of-mass subtraction. [1]


