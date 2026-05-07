import scipy.io as sio
import numpy as np

data = sio.loadmat('Params/OpenLoop_Wf_sweep.mat')
wf = data['Wf_values'].flatten()
# Let's see the names of variables in the mat file
print(data.keys())
