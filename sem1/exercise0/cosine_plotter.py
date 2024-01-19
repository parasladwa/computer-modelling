"""
Simple Python code that plots the cosine function
"""

# Import relevant python modules
import math
import matplotlib.pyplot as pyplot

# Define the cosine function
def my_function(x):
    return 

# Main method
def main():
    
    # number of data points
    n_loop = int(input('enter integer :'))

    # open output file
    file_name = str(input('enter string :'))
    out_file = open("cosine.dat","w")

    # prepare data lists
    x_values = []
    y_values = []

    # obtain function values and write them to file
    for i in range(n_loop):
        x = 2*math.pi*i/n_loop - math.pi
        f = my_function(x)
    
        # append data to lists and output file
        x_values.append(x)
        y_values.append(f)
    
        out_file.write(str(x) + " " + str(f) + "\n")

    # plot result
    pyplot.plot(x_values,y_values)
    pyplot.suptitle('Plotting the cosine function')
    pyplot.xlabel('X')
    pyplot.ylabel('Cos(X)')
    pyplot.show()
    
    
    
    
    
    
    
    


# Execute main method
if __name__ == "__main__": main()