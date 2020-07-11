#! /usr/bin/env python3

import sys
import numpy as np

sigma=0.185815
if len(sys.argv)>1:
  sigma*=float(sys.argv[1])

s_min=-2.5
s_max=+2.5

beta = 1
s_space,fes=np.loadtxt('model-fes.data',usecols=(0,1),unpack=True)
#s_space=np.linspace(-3,3,1000)

prob=np.exp(-beta*fes)
zeta=np.trapz(prob,s_space)
prob/=zeta

num_umb=1+int((s_max-s_min)/sigma)
print(' num_umb = %d'%num_umb,file=sys.stderr)

centers=np.array([sigma*i+s_min for i in range(num_umb)])
def Gauss(s):
  return np.exp(-0.5*(s/sigma)**2)

deltaF=np.array([np.trapz(prob*Gauss(s_space-c),s_space) for c in centers])

for i in range(len(s_space)):
  print(s_space[i],prob[i]*np.sum(Gauss(s_space[i]-centers)/deltaF)/num_umb,prob[i])
