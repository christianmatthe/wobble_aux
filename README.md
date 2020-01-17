# wobble_aux
auxiliary code meant to work with megbedell/wobble on CARMENES data.

Python environment Setup Instructions (preliminary):
Includes a version of wobble, do not preinstall wobble.

•conda create –name freshenv

•conda activate freshenv

•conda install tensorflow=1.12.0

•git clone https://github.com/christianmatthe/wobble_aux

•navigate to wobble folder (Wobble19_03_2019/wobble/wobble)

•python setup.py develop  (wobble/build folder and wobble/interp/interp\_op.cpython-36m-c86\_64-linux-gnu.so must be removed to before new compile (this is the case in the initial clone. Removal is onnly required if reinstalling wobble))

•install packages from wobbles requirements.txt (numpy, scipy, matplotlib, (tensorflow),tqdm, h5py, astropy)

•additional  required  packages:   barycorrpy ,pandas, dill, numpy =1.17

•wobbletest to check basic functionality
