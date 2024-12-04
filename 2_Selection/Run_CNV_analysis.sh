nano get_chrom_fasta.sh
#!/bin/bash

# Job name:
#SBATCH --job-name=seqtk
#
# Project:
#SBATCH --account=nn9449k
#
# Wall time limit:
#SBATCH --time=00-01:00:00
#
# Other parameters:
# Turn on mail notification; type can be one of BEGIN, END, FAIL, REQUEUE or ALL
#SBATCH --mail-type=ALL --mail-user=vanessa.bieker@ntnu.no
#SBATCH --cpus-per-task=1 --mem=5G

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load seqtk/1.4-GCC-12.2.0
module list

chrom=$(awk -v n=$SLURM_ARRAY_TASK_ID 'NR == n {print; exit}' /cluster/projects/nn8052k/vanessa/reindeer_reference_genome_GCA_019903745.2/ncbi_dataset/data/GCA_019903745.2/GCA_019903745.2_ULRtarCaribou_2v2_genomic.chr)
echo $chrom > $chrom.list

seqtk subseq /cluster/projects/nn8052k/vanessa/reindeer_reference_genome_GCA_019903745.2/ncbi_dataset/data/GCA_019903745.2/GCA_019903745.2_ULRtarCaribou_2v2_genomic.fasta $chrom.list > $chrom.fasta

echo "all done"

############################

sbatch --array=1-35 get_chrom_fasta.sh
# Submitted batch job 12709375
# all done



nano GenerateDuplicatedWindowRecord.sh
#!/bin/bash

# Job name:
#SBATCH --job-name=winRec
#
# Project:
#SBATCH --account=nn9449k
#
# Wall time limit:
#SBATCH --time=00-24:00:00
#
# Other parameters:
# Turn on mail notification; type can be one of BEGIN, END, FAIL, REQUEUE or ALL
#SBATCH --mail-type=ALL --mail-user=vanessa.bieker@ntnu.no
#SBATCH --cpus-per-task=1 --mem=50G

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2022.10
module list

winsize=1000


chrom=$(awk -v n=$SLURM_ARRAY_TASK_ID 'NR == n {print; exit}' /cluster/projects/nn8052k/vanessa/reindeer_reference_genome_GCA_019903745.2/ncbi_dataset/data/GCA_019903745.2/GCA_019903745.2_ULRtarCaribou_2v2_genomic.chr)
echo $chrom

echo "Step 1: Split genome into short kmer sequences."
python /cluster/home/vanessab/install/CNVcaller/bin/0.1.Kmer_Generate.py $chrom.fasta $winsize $chrom.kmer.fa

echo "Step 2: Align the kmer FASTA (from step 1) to reference genome using blasr sawriter and in the conda blasr environment."
source ${EBROOTANACONDA3}/bin/activate
conda activate /cluster/projects/nn8052k/vanessa/conda/blasr_cnv
echo "Step 2A: Create the reference.fa.sa file:"
sawriter $chrom.fasta.sa $chrom.fasta
echo "Step 2B: Align the kmer FASTA to reference genome:"
blasr $chrom.kmer.fa $chrom.fasta --sa $chrom.fasta.sa --out $chrom.kmer.aln -m 5 --noSplitSubreads --minMatch 15 --maxMatch 20 --advanceHalf --advanceExactMatches 10 --fastMaxInterval --fastSDP --aggressiveIntervalCut --bestn 10
conda deactivate

echo "Step 3: Generate duplicated window record file."
python /cluster/projects/nn8052k/0.2.Kmer_Link_VB.py $chrom.kmer.aln $winsize $chrom.window"$winsize".link

echo "all done"

############################

sbatch --array=1-35 --mem=5G --time=00-03:00:00 GenerateDuplicatedWindowRecord.sh
# Submitted batch job 12709422
# all done


# concatenate the output files in the correct order
for i in {1..35}
do
  echo $i
  chrom=$(awk -v n=$i 'NR == n {print; exit}' /cluster/projects/nn8052k/vanessa/reindeer_reference_genome_GCA_019903745.2/ncbi_dataset/data/GCA_019903745.2/GCA_019903745.2_ULRtarCaribou_2v2_genomic.chr)
  cat $chrom.window1000.link >> allChroms_window1000.link
done
# done


# Index the reference genome
# this generated a file called "referenceDB.windowsize" that should be placed in the working dir for the next steps
nano IndexRef.sh
#!/bin/bash

# Job name:
#SBATCH --job-name=indexRef
#
# Project:
#SBATCH --account=nn8052k
#
# Wall time limit:
#SBATCH --time=00-24:00:00
#
# Other parameters:
# Turn on mail notification; type can be one of BEGIN, END, FAIL, REQUEUE or ALL
#SBATCH --mail-type=ALL --mail-user=vanessa.bieker@ntnu.no
#SBATCH --cpus-per-task=1 --mem=10G

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module list

perl /cluster/home/vanessab/install/CNVcaller/bin/CNVReferenceDB.pl -w 1000 /cluster/projects/nn8052k/vanessa/reindeer_reference_genome_GCA_019903745.2/ncbi_dataset/data/GCA_019903745.2/GCA_019903745.2_ULRtarCaribou_2v2_genomic.fasta

echo "all done"

#################################

sbatch IndexRef.sh
# Submitted batch job 12709410
# done


nano filter_bam.sh
#!/bin/bash

# Job name:
#SBATCH --job-name=filterBam
#
# Project:
#SBATCH --account=nn9449k
#
# Wall time limit:
#SBATCH --time=12:00:00
#
# Other parameters:
# Turn on mail notification; type can be one of BEGIN, END, FAIL, REQUEUE or ALL
#SBATCH --mail-type=ALL --mail-user=vanessa.bieker@ntnu.no
#SBATCH --cpus-per-task=1 --mem=2G


## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge

module load SAMtools/1.18-GCC-12.3.0

module list

sampleList=$1

# get the nth line from the filelist (sample summary)
sample=$(awk -v n=$SLURM_ARRAY_TASK_ID 'NR == n' $sampleList)

echo $sample
sampleID=$(basename $sample .bam)


# filter for MAPQ 30
# filter out PCR duplicates (0x400)
# filter out not primary alignment (0x100)
# filter out supplementary alignment (0x800)
samtools view -h -q 30 -F 0x100 -F 0x400 -F 0x800 -b $sample > filtered_bam_files/$sampleID.filtered.bam

echo "all done"

########################################

ls /cluster/work/users/vanessab/coldRein_nic/share/*bam > 2024-09-24_bam.list

mkdir filtered_bam_files
sbatch --array=1-53 filter_bam.sh 2024-09-24_bam.list
# Submitted batch job 12709462
# all done

nano Ind_RD.sh
#!/bin/bash

# Job name:
#SBATCH --job-name=indRD
#
# Project:
#SBATCH --account=nn8052k
#
# Wall time limit:
#SBATCH --time=00-24:00:00
#
# Other parameters:
# Turn on mail notification; type can be one of BEGIN, END, FAIL, REQUEUE or ALL
#SBATCH --mail-type=ALL --mail-user=vanessa.bieker@ntnu.no
#SBATCH --cpus-per-task=1 --mem=10G

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load SAMtools/1.18-GCC-12.3.0
module list

bam=$(awk -v n=$SLURM_ARRAY_TASK_ID 'NR == n {print; exit}' 2024-09-26_filtered_bam.list)
sampleID=$(basename $bam .merged.rmdup.merged.realn.filtered.bam)

echo $sampleID

bash /cluster/home/vanessab/install/CNVcaller/Individual.Process.sh -b $bam -h $sampleID -d allChroms_window1000.link -s CM056435.1

echo "all done"

##################

ls $PWD/filtered_bam_files/*bam > 2024-09-26_filtered_bam.list
sbatch --array=1-53 Ind_RD.sh
# Submitted batch job 12720659

# create a list of samples to include for CNVR
ls "$PWD"/RD_normalized/* > 2024-09-27_CNVR_detection_Nic.list

# make a list of samples to exclude from CNVR detection. Their copy numbers are deduced based on the CNVR boundaries defined by other samples. This option is applicable to the outgroup or the poor quality precious samples. An empty file means all individuals are included in the CNVR detection.
touch exclude.list


nano CNVR_detection.sh
#!/bin/bash

# Job name:
#SBATCH --job-name=CNVR
#
# Project:
#SBATCH --account=nn8052k
#
# Wall time limit:
#SBATCH --time=00-24:00:00
#
# Other parameters:
# Turn on mail notification; type can be one of BEGIN, END, FAIL, REQUEUE or ALL
#SBATCH --mail-type=ALL --mail-user=vanessa.bieker@ntnu.no
#SBATCH --cpus-per-task=1 --mem=10G

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load SAMtools/1.18-GCC-12.3.0
module list

# -r: minimum of Pearson’s correlation coefficient between the two adjacent non-overlapping windows
    # 0.5 for sample size (0, 30]
    # 0.4 for sample size (30, 50]
    # 0.3 for sample size (50, 100]
    # 0.2 for sample size (100, 200]
    # 0.15 for sample size (200, 500]
    # 0.1 for sample size (500,+∞)

sampleList=$1
excludeList=$2
r=$3

bash /cluster/home/vanessab/install/CNVcaller/CNV.Discovery.sh -l $sampleList -e $excludeList -f 0.1 -h 3 -r $r -p 2024-09-27_CNVR_detection_Nic_primaryCNVR -m 2024-09-27_CNVR_detection_Nic_mergeCNVR

echo "all done"

##################
sbatch CNVR_detection.sh 2024-09-27_CNVR_detection_Nic.list exclude.list 0.3
# Submitted batch job 12723271

nano CNVR_Genotype.sh
#!/bin/bash

# Job name:
#SBATCH --job-name=cnv_geno
#
# Project:
#SBATCH --account=nn8052k
#
# Wall time limit:
#SBATCH --time=00-24:00:00
#
# Other parameters:
# Turn on mail notification; type can be one of BEGIN, END, FAIL, REQUEUE or ALL
#SBATCH --mail-type=ALL --mail-user=vanessa.bieker@ntnu.no
#SBATCH --cpus-per-task=1 --mem=10G

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module list


python /cluster/home/vanessab/install/CNVcaller/Genotype.py --cnvfile 2024-09-27_CNVR_detection_Nic_mergeCNVR --outprefix 2024-09-27_reindeer_forNic_CNV

echo "all done"

####################

sbatch CNVR_Genotype.sh
# Submitted batch job 12738465


gzcat 2024-09-27_reindeer_forNic_CNV.vcf.gz | grep "^#" > 2024-09-27_reindeer_forNic_CNV_onlyAutosome.vcf
while read line
do
  gzcat 2024-09-27_reindeer_forNic_CNV.vcf.gz | grep $line >> 2024-09-27_reindeer_forNic_CNV_onlyAutosome.vcf
done < autosomalChroms
