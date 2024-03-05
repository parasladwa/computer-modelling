operates via terminal in following format (for spyder)

%run template.py <numstep> <dt> <particle file> <xyz outfile> <OPTIONAL out>

optimal dt = 0.01, note that numstep*dt will be number of days simulation will run for.



this code should work for any given system of particles, it will print orbital periods, apsides against the systems central body (it will calculate central body itself) where applicable.



particle files in zip are :
	solar_system.txt
	mini_system.txt