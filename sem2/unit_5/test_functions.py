import template
import matplotlib.pyplot as plt

#ICs taken from 23rd of May 2023

literature_dict = {
    "Mercury": {"period": 88.0, "perihelion" : 0.30749100752, "aphelion" : 0.466584180976},
    "Venus": {"period": 224.7, "perihelion" : 0.7185931154, "aphelion" : 0.727951537368},
    "Earth": {"period": 365.2, "perihelion" : 0.983302765352, "aphelion" : 1.01672570095},
    "Mars": {"period": 687.0, "perihelion" : 1.3817041577, "aphelion" : 1.66646756902},
    "Jupiter": {"period": 4331, "perihelion" : 4.95060522107, "aphelion" : 5.45729692477},
    "Saturn": {"period": 10747, "perihelion" : 9.07499547411, "aphelion" : 10.0703304963},
    "Uranus": {"period": 30589, "perihelion" : 18.2669712228, "aphelion" : 20.063119782},
    "Neptune": {"period": 59800, "perihelion" : 29.8874574722, "aphelion" : 30.4743642214},
    "Pluto": {"period": 90560, "perihelion" : 29.658176134, "aphelion" : 49.3048461384},
    "Moon": {"period": 27.32166, "perihelion" : 0.00239876409, "aphelion" : 0.002715907640255}
}   


def percentage_difference(measured, true):
    percent = 100*(true - measured)/true
    return percent


def numstep_finder(dt, years=1):
    days = years*370
    return round(days/dt)


def extract_data(data, particles, central_body):
    organised = {}
    for particle in particles:
        if (particle.label == central_body.label) or (particle.label in organised):
            continue
        for set in data:
            if particle.label not in set:
                continue
            organised[particle.label] = set[-3:]
    return organised
    
    
#raw -- earth only
def main_earth(dts = [0.01], numstep =37000, particle_file = 'mini_system.txt', outfile = 'test_out.xyz'):
    plot_this = [dts, []]
    for dt in dts:
        print(f"____________________{dt}____________________")
        numstep = numstep_finder(dt, 5)
        particles, data, energy_deviation, central_body = template.main(numstep, dt, particle_file, outfile)
        organised_data = extract_data(data, particles, central_body)
        for pair in data:
            if pair[0] == central_body.label and pair[1] == "Earth":
                true = literature_dict[pair[1]]['period']
                measured = organised_data["Earth"][2]
                delta = percentage_difference(true, measured)
                plot_this[1].append(delta)
    print(organised_data)
    return plot_this
# a = main_earth()
# plt.plot(a[0], a[1])
# plt.show()


#refine it -- general     good dts - .2, .1, 0.05, 0.01
def main_gen(dts = [1, 0.5, 0.49, 0.2, 0.1], particle_file = 'mini_system.txt', outfile= 'test_out.xyz'):
    
    uncertainties = {
        "Earth"     : {"period" :[], "perihelion":[], "aphelion":[]},
        "Moon"      : {"period" :[], "perihelion":[], "aphelion":[]},
        "Mercury"   : {"period" :[], "perihelion":[], "aphelion":[]},
    }
    
    for dt in dts:
        
        numstep = numstep_finder(dt)
        particles, data, energy_deviation, central_body = template.main(numstep, dt, particle_file, outfile)
        data_dictionary = extract_data(data, particles, central_body)
        #print(data_dictionary)
        
        for element in data_dictionary:
            
            #period
            true = literature_dict[element]['period']
            measured = data_dictionary[element][2]
            delta = percentage_difference(measured, true)
            print(f"period uncertainty - {central_body.label}, {element}   = {delta}%\n")
            uncertainties[element]['period'].append(delta)
            
            #perihelion
            true = literature_dict[element]['perihelion']
            measured = data_dictionary[element][0]
            delta = percentage_difference(measured, true)
            print(f"perihelion uncertainty - {central_body.label}, {element}   = {delta}%\n")
            uncertainties[element]['perihelion'].append(delta)

            #aphelion
            true = literature_dict[element]['aphelion']
            measured = data_dictionary[element][1]
            delta = percentage_difference(measured, true)
            print(f"aphelion uncertainty - {central_body.label}, {element}   = {delta}%\n\n\n")
            uncertainties[element]['aphelion'].append(delta)
            
        
    #print(uncertainties)    

    for planet in uncertainties:
        print('\n', planet)
        for i in uncertainties[planet]:
            print(uncertainties[planet][i]) 
            
            '''
            plt.title(f"{planet} - {i}")
            plt.scatter(dts, uncertainties[planet][i])
            plt.show()
            '''
# main_gen()





def main_self(dts = [10, 8, 6, 4, 2, 1, 0.5, 0.1, 0.01], particle_file = 'mini_system.txt', outfile= 'test_out.xyz'):
    
    uncertainties = {
        "Earth"     : {"period" :[], "perihelion":[], "aphelion":[]},
        "Moon"      : {"period" :[], "perihelion":[], "aphelion":[]},
        "Mercury"   : {"period" :[], "perihelion":[], "aphelion":[]},
    }
    
    for dt in dts:
        
        numstep = numstep_finder(dt, 5)
        particles, data, energy_deviation, central_body = template.main(numstep, dt, particle_file, outfile)
        data_dictionary = extract_data(data, particles, central_body)
        #print(data_dictionary)
        
        for element in data_dictionary:
            
            #period
            measured = data_dictionary[element][2]
            print(f"period uncertainty - {central_body.label}, {element}   = {measured}%\n")
            uncertainties[element]['period'].append(measured)
            
            #perihelion
            measured = data_dictionary[element][0]
            print(f"perihelion uncertainty - {central_body.label}, {element}   = {measured}%\n")
            uncertainties[element]['perihelion'].append(measured)

            #aphelion
            measured = data_dictionary[element][1]
            print(f"aphelion uncertainty - {central_body.label}, {element}   = {measured}%\n\n\n")
            uncertainties[element]['aphelion'].append(measured)
            
        
    #print(uncertainties)    

    for planet in uncertainties:
        print('\n', planet)
        
        if planet == "Earth":
            for i in uncertainties[planet]:
                print(uncertainties[planet][i]) 
                
                
                plt.title(f"{planet} - {i}")
                plt.scatter(dts, uncertainties[planet][i])
                plt.show()


main_self()