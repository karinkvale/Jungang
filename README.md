Instructions for using the model
1) Download all included directories
2) UVic_ESCM is the base UVic ESCM code. This probably should not be modified. It can be placed in a "safe" spot in your file system.
3) UVOK_1.1 is modified base code that includes microplastics. It should also not be modified and can be placed in a "safe" spot next to UVic_ESCM.
4) spinup_1851 is an example working directory for actually running the model. This directory can be placed in a "working" location in your file system.
5) Please unzip A_sulphod.nc.zip in spinup_1851/data/ to A_sulphod.nc (it was too big for github without compression but needs to be uncompressed to work)
6) To make a new working directory, copy spinup_1851 to a new directory name and modify mk.in and control.in files as required. The 'updates' directory within the spinup_1851 directory are the most recent changes to the UVOK_1.1 code.
7) All code changes should be stored in the 'updates' directory within the working directory.
8) For the working directory to be able to find the model directory, a path will need to be set for your computer environment (in .profile or .bash_profile or similar)
   Add the lines
   export PATH=$PATH:$HOME/UVic_ESCM/
   ulimit -s hard
9) Enter the working directory, compile and run the model with
   mk q
