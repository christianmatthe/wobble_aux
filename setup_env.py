#Automatically sets up a Conda environment that works with wobble usinng conda and os
import os
#potential problem: conda requires responses
#######

env_name = "test_env" # ask for this in the beginning
#######
#TODO Use requirements.txt or environment.yml for smoother setup
#os.system("python3 optimize_chunk_3.1.2.py")

os.system("conda create --name {0}".format(env_name))
os.system("conda install --name {0} os".format(env_name))
os.system("conda activate {0}".format(env_name))

#Install wobble
os.system("conda install tensorflow")
os.system("python3 wobble_19_03_2019/wobble/setup.py develop")
#Install additional required packages
os.system("conda install numpy")
os.system("conda install scipy")
os.system("conda install matplotlib")
os.system("conda install tqdm")
os.system("conda install h5py")
os.system("conda install astropy")

