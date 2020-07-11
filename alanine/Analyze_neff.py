#! /usr/bin/env python3

# Used for Fig.S2b
# Calculates the effective sample size over a range of temperatures

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse


#parser
parser = argparse.ArgumentParser(description='calculate Neff over a range of temperatures')
parser.add_argument('--temp',dest='temp',type=float,default=300,required=False,help='the simulation temperature')
parser.add_argument('--mintemp',dest='mintemp',type=float,default=300,required=False,help='the minimum temperature')
parser.add_argument('--maxtemp',dest='maxtemp',type=float,default=1000,required=False,help='the maximum temperature')
parser.add_argument('--tran',dest='tran',type=int,default=0,required=False,help='transient to be skipped')
parser.add_argument('--rep',dest='wk',type=str,default='',required=False,help='replica number')
parser.add_argument('--bck',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')

args = parser.parse_args()
temp=args.temp
min_temp=args.mintemp
max_temp=args.maxtemp
wk=args.wk
if len(wk)>0:
  wk='.'+wk
outfilename='Neff'+wk+'.data'
tran=args.tran
if tran:
  outfilename='tran'+str(tran)+'-'+outfilename
  print('  tran =',tran)
bck=args.bck
if bck:
  print('  backup: '+bck)

ene_col=3
bias_col=4
#ene,bias=np.loadtxt(bck+'Colvar.data',usecols=(1,3),unpack=True,skiprows=transient)
data=pd.read_csv(bck+'Colvar'+wk+'.data',sep='\s+',comment='#',header=None,usecols=[ene_col,bias_col],skiprows=tran)
ene=np.array(data.iloc[:,0])
bias=np.array(data.iloc[:,1])
del data

kB=0.0083144621 #kj/mol
beta=1/(kB*temp)
beta_range=np.linspace(1/(kB*min_temp),1/(kB*max_temp),300)
def weights(b):
  return np.exp(-beta*((b/beta-1)*(ene-np.mean(ene))-bias))
neff=np.array([np.sum(weights(b))**2/np.sum(weights(b)**2) for b in beta_range])

cmd=subprocess.Popen('bck.meup.sh -i '+outfilename,shell=True)
cmd.wait()
np.savetxt(outfilename,np.c_[1/(kB*beta_range),neff/len(ene)],header='temp  Neff/N')
