#Neutrality tests
module load bioinfo-tools
module load ANGSD/0.933

ANC='Ancestral_cervid_genomes_mapped_CaribouV2.fasta'
REF='GCA_019903745.2_ULRtarCaribou_2v2_genomic.fa'

#1. calculate per pop saf 
angsd -b svalbard -ref $REF -anc $ANC -minMapQ 30 -minQ 30 -out out_files_II/SVA -dosaf 1 -gl 1 -remove_bads 1 -P 10 -uniqueOnly 1 -baq 1 -C 50 -minInd 8
angsd -b Nov_zem -ref $REF -anc $ANC -minMapQ 30 -minQ 30 -out out_files_II/NOZ -dosaf 1 -gl 1 -remove_bads 1 -P 10 -uniqueOnly 1 -baq 1 -C 50 -minInd 6
angsd -b russia -ref $REF -anc $ANC -minMapQ 30 -minQ 30 -out out_files_II/RU -dosaf 1 -gl 1 -remove_bads 1 -P 10 -uniqueOnly 1 -baq 1 -C 50 -minInd 6

#2. estimate sfs
realSFS ../out_files/SVA.saf.idx -P 10 > SVA.sfs
realSFS ../out_files/RU.saf.idx -P 10 > RU.sfs
realSFS ../out_files/NOZ.saf.idx -P 10 > NOZ.sfs

#3. calculate theta per site
realSFS saf2theta ../out_files/SVA.saf.idx -sfs SVA.sfs -outname SVA
realSFS saf2theta ../out_files/RU.saf.idx -sfs RU.sfs -outname RU
realSFS saf2theta ../out_files/NOZ.saf.idx -sfs NOZ.sfs -outname NOZ

#4. Estimate for every Chromosome/scaffold
thetaStat do_stat SVA.thetas.idx -outnames SVA.theta.gz
thetaStat do_stat RU.thetas.idx -outnames RU.theta.gz
thetaStat do_stat NOZ.thetas.idx -outnames NOZ.theta.gz

#5. window estimation
thetaStat do_stat SVA.thetas.idx -win 50000 -step 10000 -outnames SVA.theta.thetasWindow.gz
thetaStat do_stat RU.thetas.idx -win 50000 -step 10000 -outnames RU.theta.thetasWindow.gz
thetaStat do_stat NOZ.thetas.idx -win 50000 -step 10000 -outnames NOZ.theta.thetasWindow.gz
