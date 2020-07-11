#! /usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse


#parser
parser = argparse.ArgumentParser(description='calculate Neff over a range of temperatures and pressures')
parser.add_argument('--temp',dest='temp',type=float,default=400,required=False,help='the simulation temperature')
parser.add_argument('--mintemp',dest='mintemp',type=float,default=350,required=False,help='the minimum temperature')
parser.add_argument('--maxtemp',dest='maxtemp',type=float,default=450,required=False,help='the maximum temperature')
parser.add_argument('--pres',dest='pres',type=float,default=5000,required=False,help='the simulation preserature')
parser.add_argument('--minpres',dest='minpres',type=float,default=0,required=False,help='the minimum preserature')
parser.add_argument('--maxpres',dest='maxpres',type=float,default=10000,required=False,help='the maximum preserature')
parser.add_argument('--nbins',dest='nbins',type=int,default=100,required=False,help='number of bins')
parser.add_argument('--tran',dest='tran',type=int,default=0,required=False,help='transient to be skipped')
parser.add_argument('--bck',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-f',dest='filename',type=str,default='Colvar.data',required=False,help='input file name')
parser.add_argument('-o',dest='outfilename',type=str,default='Neff-2D.data',required=False,help='output file name')

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
bck=args.bck
if bck:
  print('  backup: '+bck)
filename=args.filename

ene_col=1
vol_col=2
bias_col=5
data=pd.read_csv(bck+filename,sep='\s+',comment='#',header=None,usecols=[ene_col,vol_col,bias_col],skiprows=tran,dtype=np.float128)
ene=np.array(data.iloc[:,0])
vol=np.array(data.iloc[:,1])
bias=np.array(data.iloc[:,2])
del data
ene-=np.mean(ene)
vol-=np.mean(vol)

kB=0.0083144621 #kj/mol
beta=1/(kB*temp)
beta_range=np.linspace(1/(kB*min_temp),1/(kB*max_temp),nbins)
pres_range=np.linspace(min_pres,max_pres,nbins)
b,p=np.meshgrid(beta_range,pres_range)
def weights(_b,_p):
  return np.exp((beta-_b)*ene+(beta*pres-_b*_p)*vol+beta*bias)

cmd=subprocess.Popen('bck.meup.sh -i '+outfilename,shell=True)
cmd.wait()
outfile=open(outfilename,'w')
print('#beta  pres  Neff/N  #N=%d'%len(ene),file=outfile)
for i in range(nbins):
  print('    working... {:.0%}'.format(i/nbins),end='\r')
  for j in range(nbins):
    print(1/(kB*b[i,j]),p[i,j]/from_bar,np.sum(weights(b[i,j],p[i,j]))**2/np.sum(weights(b[i,j],p[i,j])**2)/len(ene),file=outfile)
  print('',file=outfile)
