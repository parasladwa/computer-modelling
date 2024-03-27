#OMUAMUA
import time
import template
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

#start date of simulation
START_DATE = datetime(2023, 5, 23)

def extract_data(data, particles, central_body):
    """given data from simulation, organises into a dictionary

    Args:
        data (list): list containing particle data; periods and apsides
        particles (list): list of objects from particle class
        central_body (particle3d object): object which is central body in simulation

    Returns:
        organised(dict): organised information with nested dictionaries
    """
    organised = {}
    for particle in particles:
        if (particle.label == central_body.label) or (particle.label in organised):
            continue
        for set in data:
            if particle.label not in set:
                continue
            organised[particle.label] = set[-3:]
    return organised



def get_positions(particle_label, positions):
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



def numstep_finder(dt, years=1):
    """finds appropriate numstep according to dt

    Args:
        dt (float): timestep
        years (int, optional): number of years to find numstep for. Defaults to 1.

    Returns:
        numstep(int): number of iterations to run simulation
    """
    days = years*370
    return round(days/dt)



def main(dts = [0.3815], particle_file = 'system_omuamua.txt', outfile= 'test_out.xyz'):
    """runs the simulation and gathers the data, then looks for oumouamouas perihelion
    stores the iteration step of which it occured and runs the simulation again which 
    stores necessary data at the timestep as needed, then perihelion, velocity, and
    neptunes orbital information is analysed.

    Args:
        dts (list, optional): timestep to run . Defaults to [0.3815].
        particle_file (str, optional): filename of particle initial data. Defaults to 'system_omuamua.txt'.
        outfile (str, optional): outfile to write positions to. Defaults to 'test_out.xyz'.
    """
    dt = dts[0]

    #calculates numsteps given for a year
    numstep = numstep_finder(dt, 10)
    particles, data, energy_deviation, central_body, positions, omuamua_information = template.main(numstep, dts[0], particle_file, outfile)
    data_dictionary = extract_data(data, particles, central_body)
    
    omuamua_perihelion = data_dictionary["'Omuamua"][1]
    
    #dictionary of particles
    global particle_dict
    particle_dict = {}
    for i, particle in enumerate(particles):
        particle_dict[particle.label] = i
        particle_dict[i] = particle.label
    
    #stores positions of each needed body
    omuamua_posns = get_positions("'Omuamua", positions)
    neptune_posns = get_positions("Neptune", positions)
    sun_posns = get_positions("Sun", positions)

 
    #looks for which timestep the perihelion occured
    for i in range(0, len(omuamua_posns[0])):
        
        omuamua = np.array([omuamua_posns[0][i], omuamua_posns[1][i], omuamua_posns[2][i]])
        sun = np.array([sun_posns[0][i], sun_posns[1][i], sun_posns[2][i]])
        distance = np.linalg.norm(sun-omuamua)
        if distance == omuamua_perihelion:
            print(f"\n\n'Omuamua was at its perihelion on {(START_DATE-timedelta(days=i*dt))}") 
       
    #runs simulation again with calculated timestep of oumuamua perihelion
    particles, data, energy_deviation, central_body, positions, omuamua_information = template.main(numstep, dts[0], particle_file, outfile, extra_out=[True, i])

    #average neptune radius of orbit
    neptune_radius = (data_dictionary['Neptune'][1] + data_dictionary['Neptune'][2])/2
    
    omuamua_positions = np.transpose(omuamua_posns)
    sun_positions = np.transpose(sun_posns)
    
    outside = True
    i=0
    #once omuamua is within neptunes average orbit
    while outside:
        distance = abs(np.linalg.norm((omuamua_positions[i][:2] - sun_positions[i][:2])))
        i+=1
        if distance < neptune_radius:
            outside = False
            i-=1

    #finds suns mass
    for p in particles:
        if p.label == "Sun":
            sun_mass = p.mass
            
    #gravitational constant
    G =  8.887692593*(10**(-10))
    escape_velocity = ((2*G*sun_mass)/distance)**(1/2)

    #days behind startdate of oumuamua crosses neptunes average orbit
    entry = i*dt
    date_of_entry = (START_DATE-timedelta(days=entry))
    #prints out all necessary information
    print(f"\n\nentry date: {i} {date_of_entry}")
    print(f"'Omuamua velocity at perihelion = {np.linalg.norm(np.array(omuamua_information['velocity']))}")
    print(f"'Omuamua escape velocity : {escape_velocity} AU / day")
    print(f"'Omuamua perihelion = {omuamua_perihelion} AU")



#runs main if this file is run directly
if __name__ == "__main__":
    main()