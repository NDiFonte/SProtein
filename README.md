# SProtein

Sprotein.py is the main file to be run, before running check file intestation for run parameters
all the md_*.py  files can be ignored, as they contain all the function employed in the script

Once Sprotein.py is done, you can look for neighbor protein residues with find_neighbors.py
If a neighbor output file already exists with the same name, the script won't overwrite it

Once the neighbors are found, you can visualize the most popular encounters with plot_absolute_frequencies.py

You can also visualize the paths with VMD generating first a VMD readable trajectory with generate_VMD_traj.py (specify in the script the number of frames of the trajectory and if frame_zero should be generated)

