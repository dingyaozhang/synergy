#!/bin/bash
#SBATCH --job-name=mat4
#SBATCH --output=mat4.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50000
#SBATCH --cpus-per-task=1
#SBATCH --time=10-24:00:00

module load MATLAB
module load R


output=$1

foo () {
  local name1=$1
  local name2=$2
  res1=`python3 scripts/seesynconpvalue.py result/brcamat-log-ppi.txt  result/brcagenes-mutppi.txt ${name1} ${name2}`
  perl scripts/divsep_4ppi.pl -i ${name1} -I ${name2} -a result/expdataname.txt
  matlab -nodesktop -nosplash -nodisplay -r "addpath('scripts/');synergyGC('cache/divide/${name1}.txt', 'result/divide/${name1}.txt');exit" > /dev/null
  matlab -nodesktop -nosplash -nodisplay -r "addpath('scripts/');synergyGC('cache/divide/${name2}.txt', 'result/divide/${name2}.txt');exit" > /dev/null
  res2=`python3 scripts/seesyndrvpvalue.py result/brcadrv-ppi.txt  result/brcadrvname-ppi.txt ${name1} ${name2}`
  resall=`Rscript scripts/ppi/GCanapvalue.R result/divide/${name1}.txt result/divide/${name2}.txt $res1 $res2`
  echo $resall
}
myps=`foo ENSP00000263967 ENSP00000428056`
myps="$myps "`foo ENSP00000269305 ENSP00000263967`
myps="$myps "`foo ENSP00000269305 ENSP00000267163`
echo $myps
myps=`echo -n $myps | sed 's/ /_/g'`
Rscript scripts/fishermethod.R $myps > $output
