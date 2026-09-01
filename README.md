Simulation files for "Far-Field Thermodynamic Drag Decomposition of Canonical Bluff Bodies and Automotive Flows" Thesis by Gery Gavindra
in the scope of the Force Decomposition Project at Cranfield University in collaboration with Ford Motor Company.

How to use these files:
1. In each file case, there are several files, which are:
    - #config.txt -> this contains the experimentation history (especially on the 2D cases)
    - simSetup_[version].java -> Java script to setup the simulation file on STAR-CCM+
    - result folder -> the raw data for the last successful simulation
    - included in the result folder is a Python script used to tidy up the raw exported data from STAR-CCM+
      This only reorders the columns and rounds up the physical time in the first column (sometimes STAR-CCM results save the time with .99999, e.g., time at 5.0s saved as 4.99999999)
    - the .sim file of the latest simulation
2. The formulations used for the simulations are included in the associated Java scripts of each case. Use this to replicate the same simulation.
   The simulation is set up so all necessary parameters are stored in the automation tab, so you can easily modify some parameters by first running the Java script.
3. The included MATLAB scripts are used to do the post-processing. This includes loading the data, calculating the Strouhal number, and plotting all the decompositions.

Currently, the result for the 3D vortex-shedding cylinder is not yet available. We will upload it as soon as it is finished, along with the Windsor body case.

Any question, feel free to contact me at ggerygavindra@gmail.com
