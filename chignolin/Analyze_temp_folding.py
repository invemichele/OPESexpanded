#! /usr/bin/env python3

# Default arguments used for the inset of Fig.3
# Calculates the folded fraction at fixed pressure for different temperatures, estimating the uncertainties with block average
# CAUTION: run Prepare_analysis.sh  before this

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse


#set columns
ene_col=1
vol_col=2
bias_col=4
basin_col=5

#parser
parser = argparse.ArgumentParser(description='reweight as a function of a CV for a given temperature and pressure')
parser.add_argument('--blocks',dest='num_blocks',type=int,default=4,required=False,help='number of blocks')
parser.add_argument('--temp',dest='temp',type=float,default=500,required=False,help='the simulation temperature')
parser.add_argument('--mintemp',dest='mintemp',type=float,default=280,required=False,help='the minimum temperature')
parser.add_argument('--maxtemp',dest='maxtemp',type=float,default=370,required=False,help='the maximum temperature')
parser.add_argument('--pres',dest='pres',type=float,default=2000,required=False,help='the simulation pressure (bar)')
parser.add_argument('--rewpres',dest='rewpres',type=float,default=1,required=False,help='the reweighting pressure (bar)')
parser.add_argument('--nbins',dest='nbins',type=int,default=50,required=False,help='number of bins')
parser.add_argument('--tran',dest='tran',type=int,default=400000,required=False,help='transient to be skipped')
parser.add_argument('--bck',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-f',dest='filename',type=str,default='all_Colvar.data',required=False,help='input file name')
parser.add_argument('-o',dest='outfilename',type=str,default='temp_folded.data',required=False,help='output file name')

args = parser.parse_args()
temp=args.temp
mintemp=args.mintemp
maxtemp=args.maxtemp
from_bar=0.06022140857
pres=args.pres*from_bar
rewpres=args.rewpres*from_bar
nbins=args.nbins
outfilename=args.outfilename
tran=args.tran
if tran:
  outfilename='tran'+str(tran)+'-'+outfilename
  print('  tran =',tran,file=sys.stderr)
  tran+=1
bck=args.bck
if bck:
  print('  backup: '+bck,file=sys.stderr)
filename=args.filename

data=pd.read_csv(bck+filename,sep='\s+',comment='#',header=None,usecols=[ene_col,vol_col,bias_col,basin_col],skiprows=tran,dtype=np.float128)
ene=np.array(data.iloc[:,0])
vol=np.array(data.iloc[:,1])
bias=np.array(data.iloc[:,2])
basin=np.array(data.iloc[:,3])
del data
ene-=np.mean(ene)
vol-=np.mean(vol)


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
  log_w=((beta-b)*ene+(beta*pres-b*rewpres)*vol+beta*bias)
  return np.exp(log_w-np.amax(log_w))
folded_fraction=np.zeros(nbins)
error=np.zeros(nbins)
blocks_neff=np.zeros(nbins)
for i in range(nbins):
  print('    working... {:.0%}'.format(i/nbins),end='\r',file=sys.stderr)
  block_w=np.zeros(num_blocks)
  folded=np.zeros(num_blocks)
  w=weights(temp_range[i])
  for n in range(num_blocks):
    mask=slice(skip+n*len_blocks,skip+(n+1)*len_blocks)
    block_w[n]=np.sum(w[mask])
    folded[n]=np.sum(w[mask]*(1-basin[mask])) #folded is basin=0
  blocks_neff[i]=np.sum(block_w)**2/np.sum(block_w**2)
  folded_av=np.average(folded,axis=0,weights=block_w)
  block_w_av=np.average(block_w,axis=0,weights=block_w)
  folded_fraction[i]=folded_av/block_w_av
  error[i]=np.sqrt(1/(blocks_neff[i]-1)*np.average((folded/block_w-folded_av/block_w_av)**2,axis=0,weights=block_w))


head='temp folder_fraction error blocks_neff #num_blocks=%g'%num_blocks
cmd=subprocess.Popen('bck.meup.sh -i '+outfilename,shell=True)
cmd.wait()
np.savetxt(outfilename,np.c_[temp_range,folded_fraction,error,blocks_neff],header=head,fmt='%-9g')
