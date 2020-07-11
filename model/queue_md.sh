#!/bin/bash

# Job Settings
jname=${PWD##*/}
ncore=10
max_t=1:00 #h:min
part=''
#part='-R "select[model==XeonGold_6150]"'
singleton="-d singleton"
#logfile=log.md
outfile=log.out

#to run locally
host=$HOSTNAME
[ $# -eq 1 ] && host=$1

# Commands
mpi_cmd="plumed ves_md_linearexpansion md_input.dat"
extra_cmd="./analyze_all.sh"

# Prepare Submission
bck.meup.sh -i $outfile
#bck.meup.sh -v $logfile |& tee -a $outfile
### if euler ###
if [ ${host:0:3} == "eu-" ]
then
  cmd="mpirun ${mpi_cmd}"
  if [ ! -z "$extra_cmd" ]
  then
    cmd="${cmd}; bsub -w \"done(${jname})\" -J after$jname -o $outfile $extra_cmd"
  fi
  submit="bsub -o $outfile -J $jname -n $ncore -W $max_t $part $cmd"
  echo -e " euler submission:\n$submit" |tee -a $outfile
### if daint ###
elif [ ${host:0:5} == "daint" ]
then
  hypt="--hint=nomultithread" #avoid hyperthreading
  sb=_sbatchme.sh
  echo -e "#!/bin/bash\nsrun $hypt ${mpi_cmd}\n${extra_cmd}" > $sb
  submit="sbatch -C mc -o $outfile -J $jname -n $ncore -t $max_t:00 $part $singleton $sb"
  echo -e " daint submission:\n$submit\n $sb:" |tee -a $outfile
  cat $sb |tee -a $outfile
### if workstation ###
else
  if [ $ncore -gt 8 ]
  then
    ncore=8
  fi
  submit="time mpirun -np $ncore ${mpi_cmd}"
  echo -e " workstation submission:\n$submit\n$extra_cmd" |tee -a $outfile
  eval "$submit &>> $outfile"
  submit="$extra_cmd" # &>> $outfile"
fi

# Actual Submission
eval $submit
