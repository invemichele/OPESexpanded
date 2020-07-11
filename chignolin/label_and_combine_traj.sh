#! /bin/bash

# combine all the walkers Colvar.data in a single file, and label it as folded (0) or unfolded (1)

end=500000
if [ $# -eq 1 ]
then
  end=$1
fi

#same criteria as https://dx.plos.org/10.1371/journal.pone.0032131
low_lim=0.1
up_lim=0.4

echo "  end = $end"
bck.meup.sh -i all_Colvar.data
echo `head -1 Colvar.0.data`" basin" > all_Colvar.data
for i in `seq 0 39`
do 
  f=Colvar.${i}.data
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
        print $0,a
      }
    }
  }' $f >> all_Colvar.data
done
echo -en "  Sorting...        \r"
sort -gsk1 all_Colvar.data > tmp.all_Colvar.data
mv tmp.all_Colvar.data all_Colvar.data
