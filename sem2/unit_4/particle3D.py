"""
CompMod Ex2: Particle3D, a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
velocity and position are [3] arrays

Author: Paras Ladwa
Number: s2188899


format ; p = Particle3D ( " A " , 1.0 , [0.0 , 0.0 , 0.0] , [0.0 , 0.0 ,
0.0])


"""




import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

    Attributes
    ----------
    label: name of the particle - str
    mass: mass of the particle - float
    position: position of the particle - nparray of floats
    velocity: velocity of the particle - nparray of floats

    Methods
    -------
    __init__
    __str__
    kinetic_energy: computes the kinetic energy
    momentum: computes the linear momentum
    update_position_1st: updates the position to 1st order
    update_position_2nd: updates the position to 2nd order
    update_velocity: updates the velocity

    Static Methods
    --------------
    read_file: initializes a P3D instance from a file handle
    total_kinetic_energy: computes total K.E. of a list of particles
    com_velocity: computes centre-of-mass velocity of a list of particles
    """
    
    

    
    
    

    def __init__(self, label, mass, position, velocity):
        """
        Initialises a particle in 3D space.

        Parameters
        ----------
        label: str
            name of the particle
        mass: float
            mass of the particle
        position: [3] float array
            position vector
        velocity: [3] float array
            velocity vector
            
            tried using abs(float()) but didnt work-- for mass
        """
        ...
        
        
        #initiallising particle properties
        self.label = str(label)
        self.mass = float(mass)
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        
        
        
    

    def __str__(self):
        """
        Return an XYZ-format string. The format is
        label    x  y  z

        Returns
        -------
        str
        """
        ...
        
        #uses an fstring to return as neccessary
        return f"{self.label} {self.position[0]} {self.position[1]} {self.position[2]}"





    def kinetic_energy(self):
        """
        Returns the kinetic energy of a Particle3D instance

        Returns
        -------
        ke: float
            1/2 m v**2
        """
        ...
        
        #using temp_mass to take magnitude of velocity // repeated throughout
        temp_mass = abs(float(self.mass))
        
        
        velocity_mag = np.linalg.norm(self.velocity) #finds the magnitude of the velocity
        kinetic_e = (1/2)*(temp_mass)*(velocity_mag)**2 #applies the KE formula
        
        
        return  kinetic_e 





    def momentum(self):
        '''

        Returns the momentum of a Particle3D instance
        -------
        momentum : np array of floats
            p = mv

        '''
        
        
        temp_mass = abs(float(self.mass))
        momentum = (self.velocity)*temp_mass
        
        
        
        return momentum
    
    
    
    

    def update_position_1st(self, dt):
        """
        updates position due to only its velocity with a given time fram
        r(t + dt) = r(t) + dt · v(t)
        """
        
        self.position = self.position +  dt*self.velocity





    def update_position_2nd(self, dt, force):
        """
        r(t + dt) = r(t) + dt · v(t) + dt^2· f(t)/2m
        updates position with its velocity and a force, which is a nparray
        """
        temp_mass = abs(float(self.mass))
        self.position = self.position + np.dot(dt,self.velocity) + np.dot((dt**2), (force/(2*temp_mass)))
        
        



    def update_velocity(self, dt, force):
        #v(t + dt) = v(t) + dt · f(t)/m
        #updates velocity due to a force
        temp_mass = abs(float(self.mass))
        
        self.velocity = self.velocity + dt*(force)/temp_mass




    @staticmethod
    def read_line(line):
        """
        Creates a Particle3D instance given a line of text.

        The input line should be in the format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>

        Parameters
        ----------
        filename: str
            Readable file handle in the above format

        Returns
        -------
        p: Particle3D
        """
        
        array = line.split()
        
        #turning necessary vals into floats
        values = [float(val) for val in array[1:]]
        
        name = array[0]
        mass = values[0]
        
        pos = values[1:4]
        pos = [float(temp) for temp in pos]
    
        vel = values[4:]
        vel = [float(temp) for temp in vel]
        
        #returns instance with values read and formatted from line
        return Particle3D(name, mass, pos, vel)





    @staticmethod
    def total_kinetic_energy(particles):
        '''

        Parameters
        ----------
        particles : list
            a list of particles

        Returns
        -------
        total_kinetic : float
            net sum of kinetic energy of system from each particle

        '''
        total_kinetic = 0
        
        for particle in particles:
            #uses .kinetic_energy method to find each particles energy then adds to total
            total_kinetic += particle.kinetic_energy()
        
        
        
        return total_kinetic





    @staticmethod
    def com_velocity(particles):
        """
        Computes the CoM velocity of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        com_vel: array
            Centre-of-mass velocity
            
            
            
            the sum of each particle's momentum (mass times velocity) divided by the total mass of the system.
            
            
        """
        #total variables
        net_momentum = 0
        net_mass =0
        
        
        #iterates through list and adds momentum (from method) and mass from each particle
        for particle in particles:
            temp_mass = abs(float(particle.mass))
            net_mass+=temp_mass
            net_momentum += particle.momentum()
        
        return  net_momentum/net_mass


    










