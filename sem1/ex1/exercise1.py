"""
Code for Exercise 1

Name: Paras Ladwa
Number: s2188899
"""
import numpy as np















def task1(x, n, b):
    ''' takes an array of floating point numbers i, integer n and bool n.
    if b is true then the sum of each x element to the nth power is returned,
    f b is false, the sum of the 2nth power is returned
    '''
    result = 0
    multiple = 1
    if b ==False:
        multiple = 2    #using a variable to cover the if statement instead of having 2 loops
       
    for i in x:     # looping though the array and adding each necessary term to the result variable
        result += i**(n*multiple)
    return result
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
def task2(m, n):
    '''
    creates 3 arrays
    a 1d array of zeros of size m
    a 1d array of natural numbers of size n
    a 2d array size[m, n] of random numbers between [0,1]
    '''
    zeros = np.zeros(m)
    numbers = np.arange(1, n+1) #arange func. to create a list of the natural numbers size n
    random = np.random.random([m, n]) #by default uses range of random numbers between [0,1]
   
    return zeros, numbers, random




















def task3a(a, b, t):
    '''
    takes parameters vectors a and b, constant t
    computes and returns 2*t (a+b)
    '''
    return 2*t*(a+b)


















def task3b(x ,y):
    '''
    takes in 2 vectors and
    computs the scalar distance between them
    '''
    return np.linalg.norm(x-y) #norm func. computes magnitude















def task4a(a, b):
    '''
    takes in two vectors and tests the vector identity below
    returning the LHS and RHS seperately to be compared
    v1 × v2 = −v2 × v1
    '''
    return np.cross(a, b), np.cross(-b, a)











def task4b(a, b, c):
    '''
    very similar to part task4a
    computing both sides of the below vector identity
    which will return the LHS and RHS seperately
    v1 × (v2 + v3) = (v1 × v2) + (v1 × v3)
    '''
    return np.cross(a, (b+c)), np.cross(a, b) + np.cross(a, c)











def task4c(a, b ,c):
    '''
    very similar to part task4a
    computing both sides of the below vector identity
    which will return the LHS and RHS seperately
    v1 × (v2 × v3) = (v1 · v3)v2 − (v1 · v2)v3
    '''
    return np.cross(a, np.cross(b, c)), np.dot(a, c)*b - np.dot(a, b)*c
















def task5(x_1, x_2, m_1, m_2):
    '''
    takes in two position vectors, x_1 and x_2
    two masses m_1 and m_2
    respectively representing attributes of point masses
    computes the vector force from mass 2 acting on mass 1
    using newtons law of gravitation
    as well as the gravitational potential between the two
    and returns them individually
    '''
   
    G = 6.6743e-11
   
    r = x_2 - x_1 #r vector between x1 and x2
    r_hat = r/(np.linalg.norm(r)) #computing the unit vector
    F = (G*m_1*m_2)*(r_hat)/(np.linalg.norm(r))**2 #applying newtons law of atrtraction
    U = -(G*m_1*m_2)/(np.linalg.norm(r)) # '' gravitation
   
    return F, U #left as fully accurate numbers, not rounded as have not been advised to YET

















def task6a(n):
    '''
    takes in an integer n
    and creates a matrix of size [n, n]
    indexes through each cell
    and assigns values of i + 2j
    '''
    arr = np.zeros([n,n]) #creates an [n, n] zero matrix

    for i in range(0, n):
        for j in range(0, n): #nested iteratoin as suggested by exercise doc. to index throguh each cell
            arr[i][j] = i+2*j #changes each 0 in the array to the necessary value
    return arr















def task6b(n):
    '''
    uses task6a to create a matrix
    then computes the sum upon the axis one
    returning this as a new 1d matrix
    '''
    M = task6a(n) #recalls task6a as suggested
    new = np.sum(M, axis =1) #creates new matrix over axis =1
    return new




















