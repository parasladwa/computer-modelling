'''
project ; unit 1
    paras ladwa
    s2188899
    17/12/2023
'''
from particle3D import Particle3D
import numpy as np





def compute_separations(particles):
    '''
    Parameters
    ----------
    particles : list
        list of particle3D instances

    Returns
    -------
    seperations : numpy array
        vector seperation of each particle3D pair
    '''
    
    #set up empty array of desired dimensions
    n = len(particles)
    seperations = np.zeros([n, n, 3])
    
    
    #iterate through particle pairs
    for i in range(0, n):
        for j in range(0, n):
            
            
            #if particle indicies are identical seperation = 0
            if i==j:
                seperations[i][i] = np.zeros([3])
            
            #if pair already computed then opposite pair
            #is negative of existing vector [N3]
            elif np.linalg.norm(seperations[j][i]) != 0:
                seperations[i][j] = -seperations[j][i]
            
            
            #compute particle pair
            else:
                p1 = particles[i]
                p2 = particles[j]
                
                r = p1.position - p2.position
                
                seperations[i][j]=r
            
    return seperations
            
            



def compute_forces_potential(particles, separations):
    '''
    Parameters
    ----------
    particles : list
        list of particle3D instances
    separations : np array
        array[i][j] is seperation between
        particle[i] and particle[j] in vector form

    Returns
    -------
    forces : list
        forces[i] is net force on ith particle
        due to every other particle in particles
    
    potential_val : float
        total potential of system of particles
    '''
    n = len(particles)
    
    #gravitational constant
    G = 8.887692593*(10**(-10))
    
    
    #set up empty arrays of chosen dimensions
    potentials = np.zeros([n, n])
    forces = np.zeros([n, 3])
    
    #setting up lambda functions to reduce repetition
    mass_product = lambda i, j : particles[i].mass * particles[j].mass
    separation_norm = lambda i, j : abs(np.linalg.norm(separations[i][j]))
    
    #nested lambda functions to find forces and potentials
    potential = lambda i, j : -G * mass_product(i, j) / separation_norm(i, j)
    force = lambda i, j : potential(i, j) * separations[i][j] / (separation_norm(i, j)**2)
    
    #iterating through particle pairs
    for i in range(0, n):
        for j in range(0, n):
            
            
            #no calculations needed for particle on itself
            if i == j :
                continue
            
            #if opposite pair data is computed
            #potential is identical value[N3]
            elif potentials[j][i] != 0:
                potentials[i][j] = potentials[j][i]
                
                #compute force using lambda func
                forces[i] += force(i, j)
            
            #if opposite pair not computed
            #use lambda functions to compute
            else:
                potentials[i][j] = potential(i, j)
                
                #compute force using lambda func
                forces[i] += force(i, j)

    #sum potentials over both axis then half to account
    #for opposing pair
    potential_val = sum(sum(potentials))*0.5
    
    #return data
    return forces, potential_val
    








#made for my own testing
def make_particles():
    '''
    Returns
    -------
    particles : list of particle3D instances
        used for testing
    '''
    
    particles = [
        Particle3D("A", 1, [0, 0, 0], [0, 0, 0]),
        Particle3D("B", 1, [0, 1, 0], [0, 0, 0]),
        Particle3D("C", 1, [0, 0, 1], [0, 0, 0])
    ]
    
    return particles







