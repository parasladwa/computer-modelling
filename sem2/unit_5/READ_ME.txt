UNIT 5 READ ME FILE


AUTHOR
	paras ladwa
	s2188899



FILES
	PYTHON3 SCRIPTS
	basic_functions.py - contains smaller functions which are imported into template
	convergence.py - runs simulation repeatedly reducing dt until convergence is found
	mini_task.py - finds omuamua apsides, perihelion information
	particle3dD.py - particle class
	template.py - main simulation file

	TEXTFILES
	mini_omuamua.txt - toy system with omuamua
	mini_system.txt - toy solar system
	solar_system.txt - total solar system
	system_omuamua.txt - total solar system with omuamua

	XYZ FILE
	outfile.xyz - positions of bodies at each timestep



HOW TO

	MAIN SIMULATION

		this is within the template.py file.

		can run as is.
		parameters can be changed in main function (line 63)

		if plots are desired; bool variable in first 
		element of extra_out parameter in main (line 64).
		
		further plots are commented out if needed (lines 223 onwards).



	CONVERGENCE CALCULATIONS
	
		this is within the convergence.py file

		can run as is, it will repeatedly run simulation (10yrs) with
		decreasing dts whilst reporting to the commandline
		which observables have converged or not.

		produces plots by default, of each bodies perihelion,
		perihelion and period with every dt.
		(can be changed to False in like 220 if plots are undesired)
		
		finally prints timestep of which convergence has occured.

		number of years run for each simulation can
		be changed (default = 10yrs) on line 102



	MINI_TASK - 'Oumuamua
		
		this is within the mini_task.py file

		can run as is 
			
		number of years run for each simulation can
		be changed (default = 10yrs) on line 86

		prints out;
			date at perihelion (in first run)
			perihelion
			escape velocity at perihelion
			date of entry within Neptune's orbit 














