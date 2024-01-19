# TEST FUNC
'''

METHOD                  STATUS

init                    success
str                     success
kinetic_energy          success
momentum                success
update_position_1st     success
update_position_2nd     success
update_velocity         success
read_line               success
total_kinetic_energy    success
com_velocity            success

'''


import numpy as np
from particle3D import Particle3D



def formatting(line):
    array = line.split()
    
    values = [float(val) for val in array[1:]]
    
    name = array[0]
    mass = values[0]
    pos = values[1:4]
    vel = values[4:]
    
    array = [name, mass, pos, vel]
    
    return Particle3D(name, mass, pos, vel), array



def pull(i):    
    file = open('p3d_test_data.txt', 'r')
    f = file.readlines()
    return f[i]
        
    
    file.close()





def main():
    
    
    
    ''' dont fucking touch this - same as mess'''
    def make_particle(index):
        particle_array=[]
        for k in range(1, 6):
            particle_array += formatting(pull(k))[1]
        data = [particle_array[i:i+4] for i in range(0, len(particle_array), 4)]
        temp_particle = Particle3D(data[index][0], data[index][1], data[index][2], data[index][3])
        return temp_particle
    ''' '''

    
    
    
    
    particle = make_particle(0)
    print(f"some testing using {particle.label}")
    
    
    '''Testing'''
    
    def attributes():
        subject = particle.kinetic_energy()
        print(subject)
        
    
    
    def update_pos_1st():
        dt = float(1)
        
        print(f"befpre\n{particle.position} posn\n{particle.velocity} vel")
        
        particle.update_position_1st(dt)
        print(f"\n\nafter\n{particle.position}")
    

    def update_pos_2nd():
    
        force = np.array([1.0, 0.0, 0.0])
        print(f"2NDUPDATED TEST - {particle.label} {particle.mass} kg\n{particle.position} posn\n{particle.velocity} vel\n{force} force")
              
        particle.update_position_2nd(1, force)
        print(f"\n{particle.position} pos new")

    
    
    
    
    
    
    def array_of_particles(num):
        particles_all = []
        for j in range(0, num):
            particles_all.append(make_particle(j))
        return particles_all
    
    
    def test_total_kinetic_energy():
        particles_all = array_of_particles(5)
        return Particle3D.total_kinetic_energy(particles_all)
    
    
    def test_com_velocity():
        particles_all = array_of_particles(5)
        return Particle3D.com_velocity(particles_all)
        
    
    def test_readline():
        return Particle3D.read_line("paras 10 1 -1 0 2 -2 0")
    
main()
    

    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    