#! /usr/bin/env python3

# Used for Fig.3
# Calculates the folded fraction as a function of temperature and pressure

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

# CAUTION: run label_and_combine_traj.sh  before this, to combine the trajectories
#set columns
ene_col=1
vol_col=2
bias_col=4
basin_col=6

#parser
parser = argparse.ArgumentParser(description='calculate fraction folded and DeltaG over a range of temperatures and pressures')
parser.add_argument('--temp',dest='temp',type=float,default=500,required=False,help='the simulation temperature')
parser.add_argument('--mintemp',dest='mintemp',type=float,default=270,required=False,help='the minimum temperature')
parser.add_argument('--maxtemp',dest='maxtemp',type=float,default=800,required=False,help='the maximum temperature')
parser.add_argument('--pres',dest='pres',type=float,default=2000,required=False,help='the simulation preserature')
parser.add_argument('--minpres',dest='minpres',type=float,default=1,required=False,help='the minimum preserature')
parser.add_argument('--maxpres',dest='maxpres',type=float,default=4000,required=False,help='the maximum preserature')
parser.add_argument('--nbins',dest='nbins',type=int,default=50,required=False,help='number of bins')
parser.add_argument('--tran',dest='tran',type=int,default=0,required=False,help='transient to be skipped')
parser.add_argument('--bck',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-f',dest='filename',type=str,default='all_Colvar.data',required=False,help='input file name')
parser.add_argument('-o',dest='outfilename',type=str,default='fraction_folded.data',required=False,help='output file name')

args = parser.parse_args()
temp=args.temp
min_temp=args.mintemp
max_temp=args.maxtemp
from_bar=0.06022140857
pres=args.pres*from_bar
min_pres=args.minpres*from_bar
max_pres=args.maxpres*from_bar
nbins=args.nbins
outfilename=args.outfilename
tran=args.tran
if tran:
  outfilename='tran'+str(tran)+'-'+outfilename
  print('  tran =',tran)
  tran+=1 #first line is a comment
bck=args.bck
if bck:
  print('  backup: '+bck)
filename=args.filename

data=pd.read_csv(bck+filename,sep='\s+',comment='#',header=None,usecols=[ene_col,vol_col,bias_col,basin_col],skiprows=tran,dtype=np.float128)
ene=np.array(data.iloc[:,0])
vol=np.array(data.iloc[:,1])
bias=np.array(data.iloc[:,2])
basin=np.array(data.iloc[:,3])
del data
ene-=np.mean(ene)
vol-=np.mean(vol)
print('  all data loaded')

# f=folded, u=unfolded
f_ene=ene[basin==0]
u_ene=ene[basin==1]
if len(ene)!=len(f_ene)+len(u_ene):
  sys.exit('basin column should contain only 0 or 1')
del ene
f_vol=vol[basin==0]
u_vol=vol[basin==1]
del vol
f_bias=bias[basin==0]
u_bias=bias[basin==1]
del bias

kB=0.0083144621 #kj/mol
beta=1/(kB*temp)
temp_range=np.linspace(min_temp,max_temp,nbins)
pres_range=np.linspace(min_pres,max_pres,nbins)
t,p=np.meshgrid(temp_range,pres_range)
def f_weights(_b,_p):
  return np.exp((beta-_b)*f_ene+(beta*pres-_b*_p)*f_vol+beta*f_bias)
def u_weights(_b,_p):
  return np.exp((beta-_b)*u_ene+(beta*pres-_b*_p)*u_vol+beta*u_bias)

cmd=subprocess.Popen('bck.meup.sh -i '+outfilename,shell=True)
cmd.wait()
outfile=open(outfilename,'w')
print('#temp  pres  fraction_folded  deltaG  # N_fold=%d N_unfold=%d'%(len(f_ene),len(u_ene)),file=outfile)
for i in range(nbins):
  for j in range(nbins):
    print('    working... {:.0%}'.format((i*nbins+j)/(nbins*nbins)),end='\r')
    f_count=np.sum(f_weights(1/(kB*t[i,j]),p[i,j]))
    u_count=np.sum(u_weights(1/(kB*t[i,j]),p[i,j]))
    print(t[i,j],p[i,j]/from_bar,f_count/(f_count+u_count),-np.log(u_count/f_count),file=outfile)
  print('',file=outfile)
