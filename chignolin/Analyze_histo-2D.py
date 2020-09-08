#! /usr/bin/env python3

# Default arguments used for Fig.2a
# Calculates the histogram of the sampled target distribution as a function of energy and volume
# CAUTION: run Prepare_analysis.sh  before this

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

#parser
parser = argparse.ArgumentParser(description='calculate Histogram of energies and volumes')
parser.add_argument('--nbins',dest='nbins',type=int,default=300,required=False,help='number of bins')
parser.add_argument('--tran',dest='tran',type=int,default=400000,required=False,help='transient to be skipped')
parser.add_argument('--bck',dest='bck',type=str,default='',required=False,help='backup prefix, e.g. \"bck.0.\"')
parser.add_argument('-f',dest='filename',type=str,default='all_Colvar.data',required=False,help='input file name')
parser.add_argument('-o',dest='outfilename',type=str,default='Histo-2D.data',required=False,help='output file name')

args = parser.parse_args()
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
data=pd.read_csv(bck+filename,sep='\s+',comment='#',header=None,usecols=[ene_col,vol_col],skiprows=tran,dtype=np.float64)
ene=np.array(data.iloc[:,0])
vol=np.array(data.iloc[:,1])
del data

histo,xedges,yedges=np.histogram2d(ene,vol,nbins)
max_histo=np.max(histo)
histo/=max_histo
xcenters = (xedges[:-1] + xedges[1:]) / 2
ycenters = (yedges[:-1] + yedges[1:]) / 2
ene_mesh,vol_mesh=np.meshgrid(xcenters,ycenters)

cmd=subprocess.Popen('bck.meup.sh -i '+outfilename,shell=True)
cmd.wait()
outfile=open(outfilename,'w')
print('#ene  vol  histo  #N=%d  max_histo=%d'%(len(ene),max_histo),file=outfile)
for i in range(nbins):
  print('    working... {:.0%}'.format(i/nbins),end='\r')
  for j in range(nbins):
    print(ene_mesh[i,j],vol_mesh[i,j],histo[j,i],file=outfile)
  print('',file=outfile)
