import sys, numpy as np
p=sys.argv[1]
a=np.load(p, allow_pickle=True)
print('shape=',a.shape,'dtype=',a.dtype)
r=a.reshape(a.shape[0], -1)
print('first=', r[0][:6])
print('last=', r[-1][:6])
