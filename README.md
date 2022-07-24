
Most of the code in this repository is independent of other files and some configurations have to be done inside the files themselves.
Prerequisites:
    1. python with obspy and geographiclib (preferably installed as a conda env)
        https://github.com/obspy/obspy/wiki#installation
        https://geographiclib.sourceforge.io/

        our implementation is tested on python 3.8.12
    2. Matlab with Parallel Computing Toolbox. (can be run without but might be painfully slow for large amounts of estimation)

*code might be compatible with octave but it is untested.
Usage for real data:
    1. run the soreq.py to parse the data (you can contant me for a link to download the data)
    2. run ML_simulation_real_white.m matlab file
    3. plot the results using prints_simulation_real.m file
You can also use the .job files to run step 2 on a Slurm cluster.
You can also run with Matlab headless mode so it doesn't interrupt your workflow, using the command from a terminal inside this directory:
matlab -nosplash -nodisplay -nodesktop -sd $(pwd) -batch "ML_simulation_real_white; quit"

Usage for Synthetic Data:
    1. run F_1vs2.m file.
    2. run prints_F_per_reg.m to print the simulation results.

To use the script download obspy, geographiclib:





