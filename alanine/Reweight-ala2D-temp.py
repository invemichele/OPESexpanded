#! /usr/bin/env python3

# Used for Fig.1c
# Calculates the 2D FES over the torsion alngles, at given temperature

import sys
import numpy as np
import pandas as pd
import subprocess

#toggles
temp=300 #reweight at this temperature
temp0=300 #simulation was at this temperature
sigma=0.15
wk=''
bck=''
#bck='bck.0.'
if len(sys.argv)>1:
  temp=float(sys.argv[1])
if len(sys.argv)>2:
  wk='.'+sys.argv[2]
  print('  replica '+wk,file=sys.stderr)
print_to_file=True
if len(sys.argv)>3:
  print_to_file=False

print(' Reweighting to T =',temp)
outfilename='FES_rew2D-T'+str(temp)+'.data'

#setup
Kb=0.0083144621 #kj/mol
kbt=Kb*temp0
beta0=1./kbt
beta=1./(Kb*temp)
grid_min=-np.pi
grid_max=np.pi
grid_bin=100
cv_grid=np.linspace(grid_min,grid_max,grid_bin)
x,y=np.meshgrid(cv_grid,cv_grid)
period=grid_max-grid_min

#get kernels
filename=bck+'Colvar'+wk+'.data'
x_col=1
y_col=2
ene_col=3
bias_col=4
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[x_col,y_col,ene_col,bias_col])
cv_x=np.array(data.iloc[:,0])
cv_y=np.array(data.iloc[:,1])
ene=np.array(data.iloc[:,2])
bias=np.array(data.iloc[:,3])
del data
ene-=np.mean(ene) #numerically more stable

#build fes
basinA=0
basinB=0
max_prob=0
prob=np.zeros((grid_bin,grid_bin))
for i in range(grid_bin):
  print('    working... {:.0%}'.format(i/grid_bin),end='\r',file=sys.stderr)
  for j in range(grid_bin):
    dx=np.absolute(x[i,j]-cv_x)
    dy=np.absolute(y[i,j]-cv_y)
    arg2=(np.minimum(dx,period-dx)/sigma)**2+(np.minimum(dy,period-dy)/sigma)**2
    prob[i,j]=np.sum(np.exp((beta0-beta)*ene+bias/kbt)*np.exp(-0.5*arg2))
    if prob[i,j]>max_prob:
      max_prob=prob[i,j]
    if x[i,j]>0:
      basinB+=prob[i,j]
    else:
      basinA+=prob[i,j]

deltaF=np.log(basinA/basinB)
basinA=0
basinB=0
#deltaF_AB can be calucated also in this way (statistically compatible)
for t in range(len(bias)):
    if cv_x[t]>0:
      basinB+=np.exp((beta0-beta)*ene[t]+bias[t]/kbt)
    else:
      basinA+=np.exp((beta0-beta)*ene[t]+bias[t]/kbt)

#print out
print('  DeltaF_AB= %g DeltaF_ABbis= %g temp= %g'%(deltaF,np.log(basinA/basinB),temp))
if print_to_file:
  print('    printing...    ',end='\r',file=sys.stderr)
  with open(outfilename,'w') as f:
    print('#cv_grid  beta*FES #DeltaF_AB= %g DeltaF_ABbis= %g temp= %g'%(deltaF,np.log(basinA/basinB),temp),file=f)
    for i in range(grid_bin):
      for j in range(grid_bin):
        print(x[i,j],y[i,j],-np.log(prob[i,j]/max_prob),file=f)
      print('',file=f)
