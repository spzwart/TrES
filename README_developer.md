# User guide installation TRES
This is a short additional document that will help guide you through the installation and compilation process of TRES in some more detail, since users may want to use TRES for slightly different purposes. For a more complete guide on how to use TRES, please have a look at the [README.md](https://github.com/amusecode/TRES/blob/main/README.md) file.

For the installation, it is important to know that TRES combines two different codes: SeBa, which is a stellar evolution code, and a dynamical code. The dynamical code is already integrated into TRES and thus doesn't need to be installed separately. SeBa however, needs to be installed through the AMUSE framework. We will explain how to do this shortly, but just make sure to satisfy the [AMUSE pre-requisites](https://amuse.readthedocs.io/en/latest/install/howto-install-AMUSE.html) and have installed the necessary python modules beforehand.

If you wish to use TRES solely for simulations, without making any changes to the stellar evolution code, have a look here: [TRES as user](#TRES-as-user).

If instead you do want to be able to make changes, have a look here: [TRES as developer](#TRES-as-developer). 

In case the you have access to a computer cluster, [Run TRES on cluster](#Run-TRES-on-cluster) gives some information on running TRES with slurm.

### TRES as user
First we need to install AMUSE along with SeBa. Since we don't need all the codes from the AMUSE framework, we can just install the minimal framework and then add SeBa:

```
pip install [--user] amuse-framework
pip install [--user] amuse-<seba>
```

For clarity, this will not actually install the SeBa code, so it is not possible to change the stellar physics.

TRES can simply be installed by cloning the github repository in the terminal:

```
git clone https://github.com/amusecode/TRES.git
```

Now you should be good to go!
 
### TRES as developer
If you are interested in applying changes to the stellar evolution code SeBa, the code should be locally installed. To do this, you will have to install AMUSE directly from the source code:

```
git clone https://github.com/amusecode/amuse.git
```

Then, SeBa can be installed with:

```
pip install -e .
make seba.code
```

Now, you should be able to navigate to a directory called "amuse/src/amuse/community/seba", which contains the complete evolution code. To compile the code, write (exluding the comments preceded by a "#"):

```
make clean
make               # Here we compile all the C files
cd src/SeBa/dstar
make               # Here we make the actual executable file called "SEBA"
```

It is very important to know that anytime the code is compiled, the source code will be downloaded again, meaning any changes will be overwritten. To prevent this from happening, we need to comment out a line in the Makefile. In line 47, put a "#" before "download src/SeBa".

TRES can simply be installed by cloning the github repository in the terminal:

```
git clone https://github.com/amusecode/TRES.git
```

### Run TRES on cluster
Running computationally expensive simulations on a computer cluster can save a lot of time. However, clusters work somewhat different than your personal computer. There are two points we'd like to draw your attention to.

For starters, make sure if the pre-required packages for AMUSE are already installed. The easiest way to do this is try installing AMUSE and check where eventual errors occur. Unfortunately, on a cluster you will most likely not have sudo rights, so you'll have to figure out a way to install the missing packages. 

Second, clusters work with slurm. To run a simulation you must create a bash file that can be sumbitted as a slurm job. Down below is an example bash script for the helios cluster with additional comments for clarification. If you are using a different cluster, change the script accordingly.

```
#!/bin/bash

# first, we need to ask for resources with sbatch. The documentation can be found here: https://slurm.schedmd.com/sbatch.html
# this example uses the array functionality, which will submit an array of similar jobs that can be run simultaneously

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 1:00:00
#SBATCH -w helios-cn004
#SBATCH --job-name=TRES
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --array=0-9
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

function clean_up {
    echo "### Running Clean_up ###"
    # all files should be deleted from the compute-node before job end
    # if there are more than one output files in the directory after the completion of one of the simulations, only delete the file
    # corresponding to that simulation. Otherwise, remove the whole directory
    num_files=$(ls *.hdf | wc -l)
    echo "$num_files"
    if [ "$num_files" -eq "1" ]
    then
        rm -rf "/hddstore/$USER"
        echo "Removed all files"
    else
        #rm -rf "/hddstore/$USER/$SLURM_ARRAY_TASK_ID"
        rm "/hddstore/$USER/$FILE_NAME" 
        echo 'Removed'"$FILE_NAME"''
    fi
    echo "Finished"; date
#    # - exit the script
    exit
}

# call "clean_up" function when this script exits, it is run even if
# SLURM cancels the job
trap 'clean_up' EXIT

##### pipeline below #####

# enter virtual environment
# this environment contains the necessary python packages 
. /home/fkummer/Amuse-env/bin/activate

# openmpi was not installed on helios and needed to be imported through a module 
module purge
module load openmpi/3.1.6


mkdir -p /hddstore/$USER                              # create a directory on the node to store your data
export OUTPUT_FOLDER="/hddstore/$USER"                # this simply creates an alias for the file/directory
cd $OUTPUT_FOLDER

export FILE_NAME='TRES_'"$SLURM_ARRAY_TASK_ID"'.hdf'  # for each sbatch array, define your output filename
export TRES="/home/fkummer/TRES"
python $TRES/TPS.py -n 10 --M_max 100 --M_min 15 --Q_min 0.1 --A_distr 5 --q_min 0.1 --E_max 0.9 --E_distr 3 --e_max 0.9 --e_distr 3 -z 0.0001 -f $FILE_NAME

export FOLDER_NAME="test_TRES/cpu_check"
mkdir -p /zfs/helios/filer0/$USER/$FOLDER_NAME        # create a directory where the data should be stored (for helios this is the /zfs/helios/filer0/ directory)                             
cp $FILE_NAME /zfs/helios/filer0/$USER/$FOLDER_NAME/  # copy the data from the node to the storage directory
```
