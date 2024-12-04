
# annotate the vcf file (needed for filtering sites in LD after running plink)
module load tabix/0.2.6-GCCcore-11.3.0
VCF_in=GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup.vcf # extracts input name from command line, when calling script
VCF_out=$(basename $VCF_in .vcf)_annot.vcf # output name with _annot appended
# make annotation file with this oneliner
echo preparing annotation file
cat \
<(grep '#CHROM' $VCF_in | cut -f 1-3) \
<(paste \
<(grep -v '#' $VCF_in | cut -f 1,2) \
<(grep -v '#' $VCF_in | cut -f 1,2 | sed 's/\t/_/g')) > annotation_file.tab

echo "compressing annotation file"
bgzip annotation_file.tab

# index
echo "indexing annotation file"
tabix -p vcf annotation_file.tab.gz

module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.14-GCC-11.2.0
module list

echo "annotating vcf file using annotation file"
bcftools annotate -c CHROM,POS,ID -a annotation_file.tab.gz $VCF_in -o $VCF_out

# run plink to identify sites in LD
module --quiet purge  # Reset the modules to the system default
module load PLINK/2.00a3.7-foss-2022a
module list

plink2 --vcf GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot.vcf --indep-pairwise 50 3 0.3 -allow-extra-chr -out GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_50_3_0.3

# filter the input vcf file to only include the SNPs not in LD
# output a new vcf file
plink2 --vcf GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot.vcf --recode vcf bgz -allow-extra-chr -out GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned --extract GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_50_3_0.3.prune.in

# need to change version number in vcf file otherwise vcftools will throw an error
zcat GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned.vcf.gz | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_versionFix.vcf

module --quiet purge  # Reset the modules to the system default
module load VCFtools/0.1.16-GCC-12.3.0
# filter out sites with high missingness across samples
vcftools --vcf GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_versionFix.vcf --max-missing 0.8 --recode --out GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8

# get the sample IDs as they are in the vcf file
bcftools query -l GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.recode.vcf > samples_vcf.list

# compress and index vcf file
module load tabix/0.2.6-GCCcore-11.3.0
bgzip GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.recode.vcf
tabix -p vcf GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.recode.vcf.gz

# get a fasta file for each chromosome and sample
module --quiet purge  # Reset the modules to the system default
module load BCFtools/1.18-GCC-12.3.0
module load SAMtools/1.18-GCC-12.3.0
module list

sample=$1
chrom=1
region=$(awk -v n=$chrom 'NR == n {print; exit}' /cluster/projects/nn8052k/vanessa/reindeer_reference_genome_GCA_019903745.2/ncbi_dataset/data/GCA_019903745.2/GCA_019903745.2_ULRtarCaribou_2v2_genomic.chr)

samtools faidx /cluster/projects/nn8052k/vanessa/reindeer_reference_genome_GCA_019903745.2/ncbi_dataset/data/GCA_019903745.2/GCA_019903745.2_ULRtarCaribou_2v2_genomic.fasta $region | bcftools consensus -a ? -H I -s $sample -M N GCA_019903745.2_ULRtarCaribou_2v2_genomic.all.merged.biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.recode.vcf.gz -o fasta_files/"$sample"_chrom"$chrom".fasta


# combine all samples into one alignment file for each chromosome
# running in screen "snp" on login3
for chrom in {1..34}
do
  c=1
  for i in fasta_files/*_chrom"$chrom".fasta
  do
    echo $i $c"/63" "chrom" $chrom
    sampleID=$(basename $i _chrom"$chrom".fasta)
    echo ">"$sampleID "chrom"$chrom >> ReindeerChrom"$chrom"_samplesCombined_SNPcall_biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.fasta
    cat $i | tail -n+2 >> ReindeerChrom"$chrom"_samplesCombined_SNPcall_biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.fasta
    let c=c+1
  done

  # remove all the "?" (positions that are not in the vcf file)
  # -i: in-place: Applies the changes to the input file (no need to save it to a new file)
  sed -i 's/?//g' ReindeerChrom"$chrom"_samplesCombined_SNPcall_biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.fasta
done

# concatenate the fasta files for each chromosome into one using Geneious

raxml-ng --model GTGTR+G --all --tree rand{20},pars{20} --msa ReindeerAllChroms_samplesCombined_SNPcall_biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8.fasta --seed 02072024 --threads 50 --prefix ReindeerAllChroms_samplesCombined_SNPcall_biallelic.fmissing0.99_62genomes_outgroup_annot_LDpruned_siteMissing0.8_GTGTR_bs100 --bs-trees 100
