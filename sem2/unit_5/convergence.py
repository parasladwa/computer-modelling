import time
import template
import matplotlib.pyplot as plt



converged_dict = {
        "Earth"     : {"period" :[365.26731], "perihelion":[0.9832901529933032], "aphelion":[1.0167364502944882]},
        "Moon"      : {"period" :[27.320670370370372], "perihelion":[0.0023840651245833153], "aphelion":[0.0027185638563898586]},
        "Mercury"   : {"period" :[87.96893414634147], "perihelion":[0.3074960022651284], "aphelion":[0.46669888213499644]},
    }



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


def within(true, measured, convergence_factor):
    """checks if measured value is in range of true value with convergence_factor

    Args:
        true (float): value from dt = 0.0001
        measured (float): measured value from larger dt
        convergence_factor (float): factor which gives a range of which
                                    measured must be within to count as converged

    Returns:
        bool: determines if measured has converged
    """
    maximum = true+true*convergence_factor
    minimum = true-true*convergence_factor
    
    if minimum < measured < maximum:
        return True
    
    return False





def main(dts = [0.1], particle_file = 'mini_system.txt', outfile= 'test_out.xyz', plots = False):
    """given parametrers, runs simulation and sorts data into dictionary

    Args:
        dts (list, optional): list of dts to run simulation
        particle_file (str, optional): text file to read particle information from Defaults to 'mini_system.txt'.
        outfile (str, optional): file to which simulation stores data into. Defaults to 'test_out.xyz'.

    Returns:
        measured_vales(dict): hashmap of measured values - in same format as converged_dictionary
    """
    
    #empty hashmap which will be added into for each dt
    measured_values = {
        "Earth"     : {"period" :[], "perihelion":[], "aphelion":[]},
        "Moon"      : {"period" :[], "perihelion":[], "aphelion":[]},
        "Mercury"   : {"period" :[], "perihelion":[], "aphelion":[]},
    }
    energies = []
    
    for dt in dts:
        #iterates through timesteps and stores data into dictionary
        
        start = time.time()
        
        
        #gathers data from simulation
        numstep = numstep_finder(dt, 10)
        particles, data, energy_deviation, central_body, positions, omuamua = template.main(numstep, dt, particle_file, outfile)
        energies.append(energy_deviation)
        data_dictionary = extract_data(data, particles, central_body)
        
        #adds data into measured_values for each simulation
        for element in data_dictionary:
            
            #period
            measured = data_dictionary[element][2]
            print(f"period uncertainty - {central_body.label}, {element}   = {measured}%\n")
            measured_values[element]['period'].append(measured)
            
            #perihelion
            measured = data_dictionary[element][0]
            print(f"perihelion uncertainty - {central_body.label}, {element}   = {measured}%\n")
            measured_values[element]['perihelion'].append(measured)

            #aphelion
            measured = data_dictionary[element][1]
            print(f"aphelion uncertainty - {central_body.label}, {element}   = {measured}%\n\n\n")
            measured_values[element]['aphelion'].append(measured)
        
        #prints time taken for simulation    
        end = time.time()
        print(f"dt = {dt}")
        print(f"time elapsed = {end-start}\n\n")  
    
    
    #energy plot
    if plots:
        print(dts, energies)
        plt.title('Energy deviation with varying dts')
        plt.xlabel('dts / days')
        plt.ylabel('Energy deviation')
        plt.scatter(dts, energies)
        plt.show()
    
    #PLOTS
    if plots:
        for planet in measured_values:
            print('\n', planet)
            
            for i in measured_values[planet]:
                print(measured_values[planet][i]) 
                
                
                plt.title(f"{planet} - {i}")
                plt.xlabel('dt /days')
                if i == 'period':
                    plt.ylabel('period /days')
                else:
                    plt.ylabel(f"{i} /AU")
                plt.scatter(dts, measured_values[planet][i])
                plt.axhline(y=converged_dict[planet][i], color='r', linestyle='--', label=f'value of convergence')
                plt.legend()
                plt.show()

    return measured_values



def determine_converged():
    """runs simulation with smaller and smaller dts
    until convergence is reached based off data found from simulation
    with timestep of 0.0001 . reduces dt by 10% each iteration unitil all
    variables are considered converged.  
    """
    
    #bool holding total systems convergence state
    success = False
    dts=[]
    dt = 4
    
    #factor for which convergence is determined ie 5%
    CONVERGENCE_FACTOR = 0.005
    NUM_POINTS = 9 #3 elements per object (3objects) (apsides and period)


    while success == False:
        
        #gather data from main
        measured = main([dt])
        converged_counter = 0
        dts.append(dt)
        
        for planet in measured:
            for i in measured[planet]:
                
                measured_val = float(measured[planet][i][0])
                true_val = float(converged_dict[planet][i][0])
                
                if within(true_val, measured_val, CONVERGENCE_FACTOR):
                    converged_counter+=1
                    print(f"{planet} {i} CONVERGED")
                    
                else:
                    print(f"{planet} {i}---------------NOT CONVERGED")
                     
        print(f"counter = {converged_counter}")   
        
        #checks state of this run
        if converged_counter == NUM_POINTS:
            success = True
            
        else:#reduce dt by 10%
            dt*=.9
        
    print(f"converged at {dt}")
    print(dts)
    return dts




#runs convergence rountine if this file is run
if __name__ == "__main__":
    dts = determine_converged()
    main(dts = dts, plots = True)
    print(f"\nconvergence found at {dts[-1]} days")