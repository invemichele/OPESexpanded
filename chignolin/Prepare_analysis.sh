#! /bin/bash

# This script prepares the file all_Colvar.data, to be used for further analysis:
# - the full energy is taken from GROMACS edr files (this might be fixed in the future, see https://github.com/plumed/plumed2/issues/567 )
#   NB: GROMACS energy is very similar to the PLUMED one, and it actually does not make a difference for our observables
# - each point of the trajectory is assigned to a basin by looking at the C_alpha-RMSD using the same criteria as https://dx.plos.org/10.1371/journal.pone.0032131
# - all walkers trajectories are combined and sorted according time

end=300000
if [ $# -eq 1 ]
then
  end=$1
fi

#same criteria as https://dx.plos.org/10.1371/journal.pone.0032131
low_lim=0.1
up_lim=0.4

echo "  end = $end"
bck.meup.sh -i all_Colvar.data
echo "#! FIELDS time full_ene vol pdb_rmsd opes.bias basin" > all_Colvar.data
for i in `seq 0 39`
do 
#get full energy, with tail corrections
  ene=energy.$i.xvg
  gmx_mpi energy -f chignolin.$i.edr -o $ene <<< "12 23 0" #the volume is actually the same...
  awk -v e=$end 'NR>26{if($1==e+1) exit;print $0}' $ene > tmp.$ene #caution, header length might change
#add basin labels
  col=Colvar.${i}.data
  echo -en "  $f\r"
  if [ $i -eq 6 ] || [ $i -eq 8 ] || [ $i -eq 9 ] #specific of my initial conditions
  then 
    st=0
  else
    st=1
  fi
  awk -v a=$st -v low=$low_lim -v up=$up_lim -v e=$end '{
    if($1==e+1) exit;
    if($1=="#!") skip=1;
    else{
      if(skip==1) skip=0;
      else{
        if(a==0 && $4>up) a=1;
        if(a==1 && $4<low) a=0;
        print $4,$5,a; #rmsd bias basin
      }
    }
  }' $col > tmp.$col
#put all together
  if [ `wc -l < tmp.$ene` != `wc -l < tmp.$col` ]
  then
    echo "something is wrong, energy file and colvar file have different lengths"
    exit
  fi
  paste tmp.$ene tmp.$col |awk '{print $0}' >> all_Colvar.data
  rm tmp.$ene tmp.$col
done
echo -en "  Sorting...        \r"
sort -gsk1 all_Colvar.data > tmp.all_Colvar.data
mv tmp.all_Colvar.data all_Colvar.data
