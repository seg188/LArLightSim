

This repsoitory holds the tools for the simulation and analysis of scintillation light production in LArTPCs, with all default settings corresponding
to the DUNE near detector.


__Running the Simulation__

Make sure your enviornment has available root, numpy, and scikit-learn. All of these packages are available in conda through the conda-forge channel.

The simulation is run with default settings by either passing passing the script do_sim.py to the python interpreter: 
( $ python do_sim.py ) 

or marking do_sim.py as an executable and running:
( $ chmod +x do_sim.py )
( $ ./do_sim.py )

__File I/O__

ROOT I/O is used

All file I/O is handled in read_data.py All output files will be written to the directory stored in "WRITE_DIR" at the top of this file. The data contained in this output files is handled in the same file by the "supervisor" class. 

To add additional plots or data to the output files, you may add additional data members to the supervisor class defined in read_data.py, and ensure that they are written to the file in the member function supervisor.write(). ROOT files can store any root or C++ class. 

The supervisor class is used because it is available everywhere within the simulation as a member of the simulator class. Within the simulator class in simulator.py, the supervisor class can be accessed as self.supervisor, and its members can be filled accordingly.

__Interpretting the Output__

The output consists of whichever plots are added to the supervisor class, as outlined above, as well as six ROOT TH2D objects (in the default single-module geometry). These objects correspond to the walls of the default geometry, with each bin content being the number of hits in that bin. 

For more complex geometries, more walls will be written to the output file. Ensure that each additional module considered in the geometry has a unique name (e.g. module_1, module_2,...) to prevent naming conflicts with ROOT.

__Defining A New Geometry__

A new geometry can be defined in the do_sim.py file using the classes defined in the geometry.py file. For an example of a geometry with 2 modules, please see the file 
"geo_example.py" in the examples folder.
