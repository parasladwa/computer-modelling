import simulate_edited
import time
import numpy as np
import matplotlib.pyplot as pyplot
#%run simulate_particle.py verlet nitrogen.dat output.dat 0.01





def run(mode, particle_file, out_file='output.dat', dt=0.01):
    '''

    Parameters
    ----------
    mode : str
        euler or verlet
    particle_file : stre=
        particle file name 
    out_file : str, optional
        file to write to. The default is 'output.dat'.
    dt : float, optional
        time parameter with simulate_edited.py iterates simulation
        with. The default is 0.01.

    Returns
    -------
    E_inac : float
        energy innacuracy
    wavenum : float
        wavenumber.

    '''
    E_inac, wavenum, energies = simulate_edited.main(mode, particle_file, out_file, dt)

    return E_inac, wavenum



def find_v_0(mode, particle):
    '''

    Parameters
    ----------
    mode : str
        type of integration [euler / verlet]
    particle : str
        file name of particle data.

    Returns
    -------
    energy_unc : float
        delta E / E_0.
    wavenum : float
        wave number.

    '''
    #initialise conditions
    output_file = 'output.dat'
    
    #smallest dt
    dt=0.00001
    
    energy_unc, wavenum = run(mode, particle, output_file, dt)
    return energy_unc, wavenum




def find_v0_loops():
    '''
    finds v0 for all cases using find_v_0 function
    
    RETURNS
        list of data computed from simulations
        via run() function
    '''
    
    #lists to iterate through
    particle_list = ['nitrogen.dat', 'oxygen.dat']
    modes = ['euler', 'verlet']
    
    #to store the data
    data = [[], [], [], []]
    
    #loop through iterations
    for mode in modes:
        for particle in particle_list:
            
            #temporarily assign variables
            energy, wavenum = find_v_0(mode, particle)
            
            #append to list
            data[0].append(mode)
            data[1].append(particle)
            data[2].append(wavenum)
            data[3].append(energy)
            
            #print values for each calculation
            print(f"\n\n==============={mode} {particle}===============")
            print(f"wavenum ;  {wavenum} /cm")
            print(f"energy_inaccuracy ;  {energy}\n")
            
    return data





def different_dts(dts):
    '''
    PARAMETER
        dt ; list of floats
            values of dts to iterate through
    
    iterating through every simulation, varying dts
    and calculating errors from initial energy
    and wavenumber innacuracies
    '''
    
    
    #from excel, averages of v_0 from both
    #simulations using dt = 0.00001
    v0_oxygen = 1522.89841100558
    v0_nitrogen = 2282.83078542171
    
    #putting them in a 2d list
    particle_data = [['oxygen.dat', v0_oxygen], ['nitrogen.dat', v0_nitrogen]]
    
    #modes
    modes = ['euler', 'verlet']

    #iterating through all combinations
    for mode in modes:
        for particle in particle_data:
            for dt in dts:
                
                #using run() function to simulate with iterating conditions
                E_inaccuracy, wavenum = run(mode, particle[0], 'output.dat', dt)
                
                #convert energy inaccuracy to a percentage
                e_percent = E_inaccuracy * 100
                
                #estimate and convert wavenumber inaccuracy as percent
                true_wavenum = particle[1]
                wavenum_percent = ((true_wavenum-wavenum) / true_wavenum) * 100
                
                #outputs data
                print("\n\n===============================================")
                print(f"                 mode : {mode}")
                print(f"             particle : {particle[0]}")
                print(f"                   dt : {dt}\n")
                print("---------------------results---------------------")
                print(f"wavenumber inaccuracy : {wavenum_percent}%")
                print(f"    energy inaccuracy : {e_percent}%")
                print(wavenum)
                



#list of dts to iterate through [logarithmic scale]
dts_range = [0.01, 0.001, 0.0001, 0.00001]


different_dts(dts_range)











