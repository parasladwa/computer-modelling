import template
import matplotlib.pyplot as plt


literature_dict = {
    "Mercury": {"period": 88.0, "perihelion" : 0.30749100752, "aphelion" : 0.466584180976},
    "Venus": {"period": 224.7, "perihelion" : 0.7185931154, "aphelion" : 0.727951537368},
    "Earth": {"period": 365.2, "perihelion" : 0.983302765352, "aphelion" : 1.01672570095},
    "Mars": {"period": 687.0, "perihelion" : 1.3817041577, "aphelion" : 1.66646756902},
    "Jupiter": {"period": 4331, "perihelion" : 4.95060522107, "aphelion" : 5.45729692477},
    "Saturn": {"period": 10747, "perihelion" : 9.07499547411, "aphelion" : 10.0703304963},
    "Uranus": {"period": 30589, "perihelion" : 18.2669712228, "aphelion" : 20.063119782},
    "Neptune": {"period": 59800, "perihelion" : 29.8874574722, "aphelion" : 30.4743642214},
    "Pluto": {"period": 90560, "perihelion" : 29.658176134, "aphelion" : 49.3048461384}
}


def percentage_difference(measured, true):
    percent = 100*(true - measured)/true
    return percent


def numstep_finder(dt, years=1):
    days = years*366
    return round(days/dt)

# def main(dts = [0.01, 0.2], numstep =37000, particle_file = 'mini_system.txt', outfile = 'test_out.xyz'):
#     plot_this = [[dts], []]
#     for dt in dts:
#         data = template.main(numstep, dt, particle_file, outfile)[1]
#         for pair in data:
#             if pair[0] == "Sun":
#                 true = literature_dict[pair[1]]['period']
#                 measured = pair[4]
#                 delta = percentage_difference(true, measured)
#                 plot_this[1].append([pair[1], delta])
#                 print(pair[0], pair[1], delta)
#     return plot_this


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
    
    

def main(dts = [1, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.01], numstep =37000, particle_file = 'mini_system.txt', outfile = 'test_out.xyz'):
    plot_this = [dts, []]
    for dt in dts:
        numstep = numstep_finder(dt, 2)
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
a = main()




plt.plot(a[0], a[1])
plt.show()