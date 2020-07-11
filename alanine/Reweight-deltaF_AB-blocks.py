#! /usr/bin/env python3

# Used for Fig.1b
# Calculates the free energy difference between the alanine metastable states, as a function of temperature
# error estimate via weighted block average

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

#set columns
cv_col=1
ene_col=3
bias_col=4

#parser
parser = argparse.ArgumentParser(description='reweight deltaF as a function of temperature')
parser.add_argument('--blocks',dest='num_blocks',type=int,required=True,help='number of blocks')
parser.add_argument('--temp',dest='temp',type=float,default=300,required=False,help='the simulation temperature')
parser.add_argument('--mintemp',dest='mintemp',type=float,default=300,required=False,help='the minimum temperature')
parser.add_argument('--maxtemp',dest='maxtemp',type=float,default=1000,required=False,help='the maximum temperature')
parser.add_argument('--nbins',dest='nbins',type=int,default=50,required=False,help='number of bins')
parser.add_argument('--tran',dest='tran',type=int,default=0,required=False,help='transient to be skipped')
parser.add_argument('--bck',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-f',dest='filename',type=str,default='Colvar.data',required=False,help='input file name')
parser.add_argument('-o',dest='outfilename',type=str,default='deltaF_AB.data',required=False,help='output file name')

args = parser.parse_args()
temp=args.temp
mintemp=args.mintemp
maxtemp=args.maxtemp
nbins=args.nbins
outfilename=args.outfilename
tran=args.tran
if tran:
  outfilename='tran'+str(tran)+'-'+outfilename
  print('  tran =',tran,file=sys.stderr)
  tran+=5
bck=args.bck
if bck:
  print('  backup: '+bck,file=sys.stderr)
filename=args.filename

data=pd.read_csv(bck+filename,sep='\s+',comment='#',header=None,usecols=[cv_col,ene_col,bias_col],skiprows=tran,dtype=np.float128)
cv=np.array(data.iloc[:,0])
ene=np.array(data.iloc[:,1])
bias=np.array(data.iloc[:,2])
del data
ene-=np.mean(ene)
basin=np.zeros(len(cv)) # zero if A, one if B
for i in range(len(cv)): 
  if cv[i]>0:
    basin[i]=1

kB=0.0083144621 #kj/mol
beta=1/(kB*temp)
temp_range=np.linspace(mintemp,maxtemp,nbins)
num_blocks=args.num_blocks
len_blocks=int(np.floor(len(ene)/num_blocks))
print(' len_blocks=',len_blocks,file=sys.stderr)
skip=len(ene)-num_blocks*len_blocks
if skip!=0:
  print(' +++ WARNING blocks mismatch: throwing away first %d lines'%skip)

def weights(T):
  b=1/(kB*T)
  log_w=((beta-b)*ene+beta*bias)
  return np.exp(log_w-np.amax(log_w))
deltaF=np.zeros(nbins)
error=np.zeros(nbins)
blocks_neff=np.zeros(nbins)
for i in range(nbins):
  print('    working... {:.0%}'.format(i/nbins),end='\r',file=sys.stderr)
  block_w=np.zeros(num_blocks)
  Z_A=np.zeros(num_blocks)
  Z_B=np.zeros(num_blocks)
  w=weights(temp_range[i])
  for n in range(num_blocks):
    mask=slice(skip+n*len_blocks,skip+(n+1)*len_blocks)
    block_w[n]=np.sum(w[mask])
    Z_A[n]=np.sum(w[mask]*(1-basin[mask])) #A is basin=0
    Z_B[n]=np.sum(w[mask]*basin[mask]) #B is basin=1
  blocks_neff[i]=np.sum(block_w)**2/np.sum(block_w**2)
  deltaF[i]=np.average(-np.log(Z_B/Z_A),axis=0,weights=block_w)
  error[i]=np.sqrt(1/(blocks_neff[i]-1)*np.average((-np.log(Z_B/Z_A)-deltaF[i])**2,axis=0,weights=block_w))


head='temp deltaF_AB error blocks_neff #num_blocks=%g'%num_blocks
cmd=subprocess.Popen('bck.meup.sh -i '+outfilename,shell=True)
cmd.wait()
np.savetxt(outfilename,np.c_[temp_range,deltaF,error,blocks_neff],header=head,fmt='%-9g')
