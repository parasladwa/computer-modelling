import template







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



print('------------------------------------------------------------------------------------------')



def percentage_difference(measured, true):
    percent = 100*(true - measured)/true
    return abs(percent)



def main(dts = [0.01], numstep =37000, particle_file = 'mini_system.txt', outfile = 'test_out.xyz'):
    for dt in dts:
        data = template.main(numstep, dt, particle_file, outfile)[1]
        for pair in data:
            if pair[0] == "Sun":
                true = literature_dict[pair[1]]['period']
                measured = pair[4]
                delta = percentage_difference(true, measured)
                print(pair[0], pair[1], delta)
main()