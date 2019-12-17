# wobble_aux
auxiliary code meant to work with megbedell/wobble on CARMENES data

Python environment Setup Instructions (WIP):

• conda create –name fresh env

• conda activate fresh env

• conda install tensorflow=1.12.0 (1.10.0 apparently what megbedell used to use,
1.12.0 is in current working env) (at least tensorflow 1.14.1 should fix numpy
deprecation warnings)

• git clone wobble (paper version) or use old wobble folder (see below)

• navigate to wobble folder (Wobble 19 03 2019) (build folder and wobble/interp/interp op.cpython36m-c86 64-linux-gnu.so must be removed to force new compile)

• python setup.py develop (Breaks with master and paper version of wobble ,
but works with old Wobble 19 03 2019 version on lx39,32, laptop (not tested
whether wobble actually runs without seg fault)

• does not currently work with wobble19032019 inside wobble aux (at least on
lx32)

• install packages from requirements.txt (numpy, scipy, matplotllib, (tensorflow),
tqdm, h5py, astropy)

• additional required packages: barycorrpy(apparently not conda installable?),
pandas, [conda install -c astropy astroquery], dill, (numpy=1.16 if deprecation
warning, not compatible with ”extrapolate”) numpy =1.17

• wobble test to check basic functionality
