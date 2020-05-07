#!/bin/bash
#SBATCH -A SNIC2020-7-19
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -o out.txt
#SBATCH -t 50:00:00
#SBATCH --mail-user=cheyu@chalmers.se
#SBATCH --mail-type=end
module load intel/2018b CMake GMP
i=41
a1=$i+0
b1=$i+0
a2=$i+1
b2=$i+1
a3=$i+2
b3=$i+2
a4=$i+3
b4=$i+3
a5=$i+4
b5=$i+4
a6=$i+5
b6=$i+5
a7=$i+6
b7=$i+6
a8=$i+7
b8=$i+7
a9=$i+8
b9=$i+8
a10=$i+9
b10=$i+9
a11=$i+10
b11=$i+10
a12=$i+11
b12=$i+11
a13=$i+12
b13=$i+12
a14=$i+13
b14=$i+13
a15=$i+14
b15=$i+14
a16=$i+15
b16=$i+15
a17=$i+16
b17=$i+16
a18=$i+17
b18=$i+17
a19=$i+18
b19=$i+18
a20=$i+19
b20=$i+19
matlab -nodesktop -singleCompThread -r "simulationBiolog($a1,$b1)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a2,$b2)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a3,$b3)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a4,$b4)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a5,$b5)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a6,$b6)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a7,$b7)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a8,$b8)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a9,$b9)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a10,$b10)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a11,$b11)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a12,$b12)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a13,$b13)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a14,$b14)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a15,$b15)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a16,$b16)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a17,$b17)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a18,$b18)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a19,$b19)" &
matlab -nodesktop -singleCompThread -r "simulationBiolog($a20,$b20)" &
wait;
