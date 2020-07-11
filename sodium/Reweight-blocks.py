#! /usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

#set columns
ene_col=1
vol_col=2
cv_col=4
bias_col=5
#set cv stuff
cv_min=0
cv_max=1
rescale_cv=250 #to have crystallyinity from 0 to 1

#parser
parser = argparse.ArgumentParser(description='reweight as a function of a CV for a given temperature and pressure')
parser.add_argument('--blocks',dest='num_blocks',type=int,required=True,help='number of blocks')
parser.add_argument('--temp',dest='temp',type=float,default=400,required=False,help='the simulation temperature')
parser.add_argument('--rewtemp',dest='rewtemp',type=float,required=True,help='the reweighting temperature')
parser.add_argument('--pres',dest='pres',type=float,default=5000,required=False,help='the simulation pressure (bar)')
parser.add_argument('--rewpres',dest='rewpres',type=float,required=True,help='the reweighting pressure (bar)')
parser.add_argument('--sigma',dest='sigma',type=float,default=0.01,required=False,help='sigma for KDE')
parser.add_argument('--nbins',dest='nbins',type=int,default=100,required=False,help='number of bins')
parser.add_argument('--tran',dest='tran',type=int,default=0,required=False,help='transient to be skipped')
parser.add_argument('--bck',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-f',dest='filename',type=str,default='all_Colvar.data',required=False,help='input file name')
parser.add_argument('-o',dest='outfilename',type=str,default='FES_rew.data',required=False,help='output file name')
parser.add_argument('--nomintozero',dest='nomintozero',action='store_true',default=False,help='do not shift the minimum to zero')

args = parser.parse_args()
temp=args.temp
rewtemp=args.rewtemp
from_bar=0.06022140857
pres=args.pres*from_bar
rewpres=args.rewpres*from_bar
sigma=args.sigma
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

data=pd.read_csv(bck+filename,sep='\s+',comment='#',header=None,usecols=[ene_col,vol_col,cv_col,bias_col],skiprows=tran,dtype=np.float128)
ene=np.array(data.iloc[:,0])
vol=np.array(data.iloc[:,1])
cv=np.array(data.iloc[:,2])
bias=np.array(data.iloc[:,3])
del data
ene-=np.mean(ene)
vol-=np.mean(vol)
cv/=rescale_cv

kB=0.0083144621 #kj/mol
beta=1/(kB*temp)
rewbeta=1/(kB*rewtemp)
cv_grid=np.linspace(cv_min,cv_max,nbins)

cmd=subprocess.Popen('bck.meup.sh -i '+outfilename,shell=True)
cmd.wait()

num_blocks=args.num_blocks
len_blocks=int(np.floor(len(cv)/num_blocks))
if num_blocks*len_blocks!=len(cv):
  print(' +++ WARNING blocks mismatch: throwing away last %d lines'%(len(cv)-num_blocks*len_blocks))

weight=np.exp((beta-rewbeta)*ene+(beta*pres-rewbeta*rewpres)*vol+beta*bias)
block_w=np.zeros(num_blocks)
prob=np.zeros((num_blocks,nbins))
for n in range(num_blocks):
  for i in range(nbins):
    print('    working... {:.0%}'.format((n*nbins+i)/nbins/num_blocks),end='\r')
    block_w[n]=np.sum(weight[n*len_blocks:(n+1)*len_blocks])
    prob[n,i]=np.sum(weight[n*len_blocks:(n+1)*len_blocks]*np.exp(-0.5*((cv_grid[i]-cv[n*len_blocks:(n+1)*len_blocks])/sigma)**2))
blocks_neff=np.sum(block_w)**2/np.sum(block_w**2)
av_fes=np.average(-np.log(prob),axis=0,weights=block_w)
blocks_var=blocks_neff/(blocks_neff-1)*np.average((-np.log(prob)-av_fes)**2,axis=0,weights=block_w)
#av_prob=np.average(prob,axis=0,weights=block_w)
#blocks_var=blocks_neff/(blocks_neff-1)*np.average((prob-av_prob)**2,axis=0,weights=block_w)
error=np.sqrt(blocks_var/blocks_neff)

#av_fes=-np.log(av_prob)
#error/=av_prob #error propagation
if not args.nomintozero:
  av_fes-=min(av_fes)

np.savetxt(outfilename,np.c_[cv_grid,av_fes,error],header='cv FES error #temp= %g K, pres= %g bar, blocks_neff=%g'%(rewtemp,rewpres/from_bar,blocks_neff))
