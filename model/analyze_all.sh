#! /bin/bash

# Used for Fig.S4

for i in `seq 0 9`
do
  echo " --replica $i"
  ./Reweight-multi.py -r $i
#  ./Reweight-multi.py -r $i -f
done


make_stats() {
  bck.meup.sh -i $2
  paste $1 |\
    awk 'BEGIN{print "#average std_dev"}
         NR>1{
           av=0; av2=0; 
           for (i=2; i<=NF; i+=2) {av+=$i; av2+=$i^2;};
           av*=2./NF; av2*=2./NF; 
           print $1,av,sqrt(av2-av^2)
         }' > $2
}

bck=''

make_stats "${bck}fes_deltaF.rew.?.data" Stats-fes_deltaF.rew.data
#make_stats "tran-1/${bck}fes_deltaF.rew.?.data" Stats-flip-fes_deltaF.rew.data
