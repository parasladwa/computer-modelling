from particle3D import Particle3D




def particles_from_file(filename = "mini_system.txt"):
    '''
    Parameters
    ----------
    filename : str, opt
        filename of particles
        . The default is "mini_system.txt".

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





particles = particles_from_file()
for p in particles:
    print(p)