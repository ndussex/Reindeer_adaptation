#step 1

#!/bin/bash -l
#SBATCH -A Proj_ID
#SBATCH -p node -n 1
#SBATCH -J PBS_Rein_step1
#SBATCH -t 7-00:00:00

module load bioinfo-tools
module load ANGSD/0.917

#snic2022-22-639
#snic2022-5-27

date
mkdir out_files/ 

REF='Ancestral_cervid_genomes.fasta'

# calculate per pop saf for each population
angsd -b svalbard -anc $REF -minMapQ 30 -minQ 30 -out out_files/SVA -dosaf 1 -gl 1
#angsd -b Nov_zem -anc $REF -minMapQ 30 -minQ 30 -out out_files/NOZ -dosaf 1 -gl 1
#angsd -b russia_II -anc $REF -minMapQ 30 -minQ 30 -out out_files/RU -dosaf 1 -gl 1


#step 2

#!/bin/bash -l
#SBATCH -A Proj_ID
#SBATCH -p node -n 1
#SBATCH -J PBS_Rein_step2
#SBATCH -t 6-00:00:00

module load bioinfo-tools
module load ANGSD/0.917

#snic2022-22-639
#snic2022-5-27

date

REF='Ancestral_cervid_genomes.fasta'

# calculate per pop saf for each population

# calculate all pairwise 2dsfs's
# grep '>' $REF | sed '/>/s///' > scf.list
# need to use the -fold 1 option if  no ancstral state

cat scf_caribouV2.list | while read LINE
do
	realSFS -P 10 -r $LINE out_files/SVA.saf.idx out_files/NOZ.saf.idx > out_files/SVA.NOZ.$LINE.ml   # PBS0
	realSFS -P 10 -r $LINE out_files/SVA.saf.idx out_files/RU.saf.idx > out_files/SVA.RU.$LINE.ml     # PBS1
	realSFS -P 10 -r $LINE out_files/RU.saf.idx out_files/NOZ.saf.idx > out_files/RU.NOZ.$LINE.ml     # PBS2
done

#date

# prepare the fst for easy analysis etc
cat scf_caribouV2.list | while read LINE
do
	realSFS fst index -r $LINE out_files/SVA.saf.idx out_files/NOZ.saf.idx out_files/RU.saf.idx -sfs out_files/SVA.NOZ.$LINE.ml -sfs out_files/SVA.RU.$LINE.ml -sfs out_files/RU.NOZ.$LINE.ml -fst
out out_files/pbs.$LINE
done

date

#print results in sliding windows
cat scf_caribouV2.list | while read LINE
do
	realSFS fst stats2 out_files/pbs.$LINE.fst.idx -win 50000 -step 10000 > out_files/pbs.$LINE.sliding_window
done

date
