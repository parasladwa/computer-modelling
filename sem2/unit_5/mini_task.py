#OMUAMUA
import time
import template
import numpy as np
import convergence
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

START_DATE = datetime(2023, 5, 23)

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



def main(dts = [0.1], particle_file = 'system_omuamua.txt', outfile= 'test_out.xyz'):
    dt = dts[0]
    start = time.time()
    
    numstep = convergence.numstep_finder(dt, 10)
    particles, data, energy_deviation, central_body, positions = template.main(numstep, dts[0], particle_file, outfile)
    data_dictionary = convergence.extract_data(data, particles, central_body)
    
    omuamua_perihelion = data_dictionary["'Omuamua"][1]
    
    print(data_dictionary)
    
    global particle_dict
    particle_dict = {}
    for i, particle in enumerate(particles):
        particle_dict[particle.label] = i
        particle_dict[i] = particle.label
    
    omuamua_posns = get_positions("'Omuamua", positions)
    neptune_posns = get_positions("Neptune", positions)
    sun_posns = get_positions("Sun", positions)

    end = time.time()
    print(f"time elapsed = {end-start}\n\n")  
    
    for i in range(0, len(omuamua_posns[0])):
        
        omuamua = np.array([omuamua_posns[0][i], omuamua_posns[1][i], omuamua_posns[2][i]])
        sun = np.array([sun_posns[0][i], sun_posns[1][i], sun_posns[2][i]])
        distance = np.linalg.norm(sun-omuamua)
        if distance == omuamua_perihelion:
            print(f"'Omuamua was at its perihelion on {(START_DATE-timedelta(days=i*dt))}")    
    
            
            
    particles, data, energy_deviation, central_body, positions = template.main(numstep, dts[0], particle_file, outfile, extra_out=[True, i])

    
        
        
    
    
    
    
main()
