# wobble_aux
auxiliary code meant to work with megbedell/wobble on CARMENES data

Python environment Setup Instructions (WIP):

•conda create –name freshenv

•conda activate freshenv

•conda install tensorflow=1.12.0 (1.10.0 apparently what megbedell used to use,1.12.0 is in current working env) (at least tensorflow 1.14.1 should fix numpydeprecation warnings)

•git clone https://github.com/christianmatthe/wobble_aux

•navigate to wobble folder (Wobble19_03_2019/wobble/wobble)

•python setup.py develop

•install packages from requirements.txt (numpy, scipy, matplotllib, (tensorflow),tqdm, h5py, astropy)

•additional  required  packages:   barycorrpy(apparently  not  conda  installable?),pandas, [conda install -c astropy astroquery], dill, (numpy=1.16 if deprecationwarning, not compatible with ”extrapolate”) numpy =1.17

•wobbletest to check basic functionality
