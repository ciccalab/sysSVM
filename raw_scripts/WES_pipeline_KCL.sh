#!/bin/sh

## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd

#$ -j y

## Job Name, can be anything##
#$ -N aln_test

## Set the SHELL type to be sh##
#$ -S /bin/sh


# -------------------------------------------------------------------------------------------------#
####################################################################################################
####################################################################################################
# -------------------------------------------------------------------------------------------------#
# START FUCTION DEFINITION

function get_dmg()
{
	
	local sample=$1
	local workDir=$2
	local outDir=$3
	local dbNSFP_path=$4
	
	mkdir -p $outDir
	
	local onT="$workDir/$sample.varscan.all.BP.OT.vcf"	
	local infile="$workDir/$sample.dmg.input.bed"
	local dbNSFP2_ver="search_dbNSFP24"
	local outFile="$outDir/$sample.dmg.bed"
	
	# make a zero based coordinate file for annovar
	awk 'NR > 1 {print $1"\t"$2"\t"$3"\t"$19}' $onT > $infile
	
	cd $dbNSFP_path
	echo "... Calling Damaging Mutations"
	java -Xmx8g $dbNSFP2_ver -i $infile -o $outFile -v hg19
	
	mv $outFile $workDir
	 
	rm $infile
	rm -r $outDir
	
	cd $workDir
}


function get_annovar()
{	

	# Annovar Config 
	local params="humandb/ -buildver hg19 -protocol refGene,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_afr,1000g2012apr_amr,1000g2012apr_asn,esp6500si_all,esp6500si_ea,esp6500si_aa,nci60,snp138,clinvar_20140211,cosmic68 -operation g,f,f,f,f,f,f,f,f,f,f,f,f -nastring NA"
	
	local annovar_cmd="perl table_annovar.pl"
		
	# Local Variables
	local annovar=$1
	local sample=$2
	local workDir=$3
	local outDir=$4
		
	local onT="$workDir/$sample.varscan.all.BP.OT.vcf"	
	local infile="$workDir/$sample.varscan.all.BP.OT.0base.bed"
	
	mkdir -p $outDir
	
	# make a zero based coordinate file for annovar
	awk 'NR > 1 {print $1"\t"$2"\t"$2"\t"$3"\t"$19}' $onT > $infile
		
	#cd $annovar
	cd "/home/FC/Software/annovar_CRC/"
	
	echo "run annovar for $sample..."
	$annovar_cmd $infile $params --outfile "$outDir/$sample"
	
	cd $workDir
	rm $infile	
}


function get_annovar_indels()
{	

	# Annovar Config 
	local params="humandb/ -buildver hg19 -protocol refGene,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_afr,1000g2012apr_amr,1000g2012apr_asn,esp6500si_all,esp6500si_ea,esp6500si_aa,nci60,snp138,clinvar_20140211,cosmic68 -operation g,f,f,f,f,f,f,f,f,f,f,f,f -nastring NA"
	
	local annovar_cmd="perl table_annovar.pl"
		
	# Local Variables
	local annovar=$1
	local sample=$2
	local workDir=$3
	local outDir=$4
		
	local onT="$workDir/$sample.varscan.indel.all.BP.OT.vcf"
	local infile="$workDir/$sample.varscan.indel.all.BP.OT.0base.bed"
	
	mkdir -p $outDir
	
	# make a the annovar file for indels 
	R CMD BATCH --vanilla "--args $workDir $sample" "$varc2annovar" "/dev/null"
		
	#cd $annovar
	cd "/home/FC/Software/annovar_CRC/"
	
	echo "run annovar for $sample..."
	$annovar_cmd $infile $params --outfile "$outDir/$sample"
	
	cd $workDir
	rm $infile	
}


function get_stat_uniCov()
{
	local fileIn=$1
	local fileOut=$2
	local s=$3
	local ftarget=$4
	
	echo " ... computing the uniformity of coverage stats for $fileIn"
	eql_zero=$( awk 'NR>1 {sum+=($3-$2)} END {print sum}' $ftarget)
	grt_zero=($(awk '$1 > 0' $fileIn | wc -l ))
	grt_10=($(awk '$1 >= 10' $fileIn | wc -l ))
	grt_30=($(awk '$1 >= 30' $fileIn | wc -l ))
	grt_60=($(awk '$1 >= 60' $fileIn | wc -l ))
	grt_100=($(awk '$1 >= 100' $fileIn | wc -l ))
	grt_300=($(awk '$1 >= 300' $fileIn | wc -l ))
	
	echo -e "sample\t=0x\t>0x\t>=10x\t>=30x\t>=60x\t>=100x\t>=300" > $fileOut
	echo -e "$s\t$eql_zero\t$grt_zero\t$grt_10\t$grt_30\t$grt_60\t$grt_100\t$grt_300" >> $fileOut
}

function get_stat_uniCov50()
{
	local fileIn=$1
	local fileOut=$2
	local s=$3
	local ftarget=$4
	
	echo " ... computing the uniformity of coverage stats for $fileIn"
	eql_zero=$( awk 'NR>1 {sum+=($3-$2)} END {print sum}' $ftarget)
	grt_zero=($(awk '$1 > 0 && $1 <= 50' $fileIn | wc -l ))
	grt_60=($(awk '$1 > 50 && $1 <= 200' $fileIn | wc -l ))
	grt_100=($(awk '$1 > 200' $fileIn | wc -l ))
	
	echo -e "sample\t=0x\t]0x,50x]\t]50x,200x]\t>200x" > $fileOut
	echo -e "$s\t$eql_zero\t$grt_zero\t$grt_60\t$grt_100" >> $fileOut
}


function count_raw_reads {
	
	printf "## Raw reads in fastq (not filtered)\n" >> $1
	IFS=$'\x0A'
	# R1
	names=( $(ls -dtr $2) )
	rawreads=0
	for file in ${names[@]}
	do
		printf "${file}\t" >> $1
		cnt=$(zcat $file | grep $4 | wc -l)
		printf "${cnt}\n" >> $1
		rawreads=$(($rawreads + $cnt))
	done;
	printf "# TOT RAW READS R1:\t ${rawreads}\n" >> $1
	printf "\n" >> $1
	
	# R2
	names=( $(ls -dtr $3) )
	rawreads=0
	for file in ${names[@]}
	do
		printf "${file}\t" >> $1
		cnt=$(zcat $file | grep $4 | wc -l)
		printf "${cnt}\n" >> $1
		rawreads=$(($rawreads + $cnt))
	done;
	printf "# TOT RAW READS R2:\t ${rawreads}\n" >> $1
	printf "\n" >> $1
}


function filer_raw_reads {
	
	local path=$1
	local ffiltered=$2
	local flog=$3
	local HISEQ=$4

	local R1_reads=( $(ls -dtr *"R1"*) )
	local R2_reads=( $(ls -dtr *"R2"*) )
	
	local old_dir=$(pwd)
	cd $path
	
	printf "## Filtered reads in fastq\n" >> $flog
	
	printf "$R1_reads\t" >> $flog
	zgrep -A 3 '^@.* [^:]*:N:[^:]*:' "$path/$R1_reads" | grep -v -- '^--$' | gzip -c > "$path/$ffiltered/filtered_$R1_reads" ;
	cnt=$(zgrep $HISEQ "$path/$ffiltered/filtered_$R1_reads" | wc -l)
	printf "${cnt}\n" >> $flog
	printf "# TOT FILTERED READS R1:\t ${cnt}\n" >> $flog
	printf "\n" >> $flog
	
	printf "$R2_reads\t" >> $flog
	zgrep -A 3 '^@.* [^:]*:N:[^:]*:' "$path/$R2_reads" | grep -v -- '^--$' | gzip -c > "$path/$ffiltered/filtered_$R2_reads" ;
	cnt=$(zgrep $HISEQ "$path/$ffiltered/filtered_$R2_reads" | wc -l)
	printf "${cnt}\n" >> $flog
	printf "# TOT FILTERED READS R2:\t ${cnt}\n" >> $flog
	printf "\n" >> $flog
	
	cd $old_dir
	echo $old_dir
}

#   END FUNCTION DEFINITION                                                                        #
# -------------------------------------------------------------------------------------------------#
####################################################################################################
####################################################################################################
# -------------------------------------------------------------------------------------------------#
#                                                                                        #


  
# CONFIGURATION
# --------------

sample=$1 # Name of the sample in input

# Libraries, modules and tools
export PATH=~/bin_tools/samtools-0.1.18:$PATH
export PATH=~/bin_tools/bedtools_2.19.1:$PATH
export PATH=~/bin_tools/R-3.0.3/bin:$PATH
module load novocraft/3.01.02
varscn='/home/gambardg/bin_tools/VarScan.v2.3.6.jar'
annovar_path='/home/FC/Software/annovar_CRC/'
dbsnfs="/home/FC/DB/dbNSFP2.4"


# where are the fastq files
path="/home/FC/ClonalExpansion/WholeExome/150129_KCL/$sample"

# Where the results will be stored
root="/home/FC/ClonalExpansion/WholeExome/150129_KCL/$sample"

# library config
mean_len_library=400
sd_len_library=100
HISEQ='HISEQ'

# On target configurations:
ftarget='/home/FC/Capture_kits/Human_V5/S04380110/S04380110_Covered.bed'

# ref genome
genome_ref='/home/FC/DB/Genomes/Hs_v37/hg19.fa'
hg19_index="/home/FC/DB/hg19/hg19.13.2.nix"

# samtools
samtls='samtools'

# R scripts
varc2annovar="/home/FC/ClonalExpansion/WholeExome/varc2annorINDELs.R"


# END configuration


# MAIN
# ------

# aln folder
align='aligned'

# filtered fastq folder
ffiltered='filtered-fastq'

OLD_IFS=$IFS
mkdir -p $root
cd $root

mkdir -p "$path/$ffiltered"
mkdir -p "$root/$align"


# 00. RAW READS
# ---------------

current=$(date +"%Y%m%d")
flog="${root}/$align/${current}_$sample.log"

# We are assunming ONLY two fastq files!!!
echo ".. Counting Reads (Step 00)"
R1_reads=( $(ls -dtr "$path/"*"R1"*) )
R2_reads=( $(ls -dtr "$path/"*"R2"*) )
count_raw_reads $flog $R1_reads $R2_reads $HISEQ

# 01. FILTERING RAW READS
# ------------------------

# We are assunming ONLY two fastq files!!!
echo ".. Filtering Reads (Step 01)"
filer_raw_reads $path $ffiltered $flog $HISEQ

# 02. ALIGNMENT
# ---------------

# we move in the aln directory
cd "$root/$align"

fout="$sample.sam"

# get the fastq filterd file paths
R1_reads=( $(ls -dtr "$path/$ffiltered/"*"R1"*) )
R2_reads=( $(ls -dtr "$path/$ffiltered/"*"R2"*) )

# Alignment with novalign
echo ".. Allinment and stats (Step 02)"
novoalign -o SAM -d $hg19_index -F STDFQ -i PE $mean_len_library,$sd_len_library -f $R1_reads $R2_reads > $fout

# Mismatches <=3 e Uniquely Aligned and Sorted
perl -e 'while (<>) {@f=split/\s+/; next if ($f[5]=~/N|H|P/); if ($f[15]=~/34/) {print;next} $a=($f[15]=~s/\d+[A-Z]//g); print if $a<=3}' $fout | $samtls view - -S -F 4 -b -u | $samtls sort - "${fout}.srt"

# Counts ALIGNED reads with Dups
$samtls flagstat $fout.srt.bam > "$sample.aligned.stat.txt"


# Remove dublicatesf
$samtls rmdup $fout.srt.bam $fout.bam

# Counts ALIGNED reads NO Dups
$samtls flagstat $fout.bam > "$sample.aligned.nodup.stat.txt"


# conto READS ONTARGET che comprendono anche le regioni overlappanti
printf "# TOT (ALIGNED) READS ON TARGET:\t " >> $flog
intersectBed -abam $fout.bam -b $ftarget > "${fout}_srt_OT.bam";
$samtls view "${fout}_srt_OT.bam" | grep $HISEQ | wc -l >> $flog
printf "\n" >> $flog

$samtls flagstat "${fout}_srt_OT.bam" > "$sample.aligned.nodup.stat.BP.OT.txt"
$samtls index "${fout}_srt_OT.bam"


printf "# TOT (ALIGNED) BASES ON TARGET:\t " >> $flog
$samtls mpileup -f $genome_ref $fout.bam | awk -F "\t" '{s=$2-1; print $1"\t"s"\t"$2"\t"$4}' | intersectBed -a stdin -b $ftarget > "${fout}_BP_OT.bed"
perl -ane '$s+=$F[3];print $s."\n"' "${fout}_BP_OT.bed" | tail -1 >> $flog #perl -ane '$s+=$F[4];print $s."\n"' "${i}_BP_OT.bed" | tail -1 >> $flog
cut -f 4 "${fout}_BP_OT.bed" > "${fout}_BP_OT.cov"
printf "\n" >> $flog

fbam="${i}_${j}_srt_OT.bam"
coverageBed -abam $fbam -b $ftarget > "${i}_${j}.covarage"
perl -ane 'print if $F[5]==0' "${i}_${j}.covarage" > targets_not_covered.bed

# Get uniformity of coverage
cd "$root/$align"
get_stat_uniCov "$sample.sam_BP_OT.cov" "$sample.sam_BP_OT.cov.txt" $sample $ftarget
get_stat_uniCov50 "$sample.sam_BP_OT.cov" "$sample.sam_BP_OT.cov.50x.txt" $sample $ftarget



# 03. Variant Calling
# --------------------

# 00a. Configuration: necessary file path
# ----------------------------------------#
root="$root/$align"

# 00b. Configuration: Varscan Parameters
# ----------------------------------------#
minCov=1
minVar=1
minQ=15


# Let be sure to be in the right directory
cd $root



# Make FIFO pileups (SNP + INDELs)
fname="$sample.pileup.fifo"
mkfifo $fname

fname_indel="$sample.pileup.indel.fifo"
mkfifo $fname_indel

fbam="$sample.sam_srt_OT.bam"
if [ -e $fbam ]
	then
		samtools mpileup -q 1 -d 10000 -f $genome_ref $fbam > $fname &
		samtools mpileup -q 1 -d 10000 -f $genome_ref $fbam > $fname_indel &
		echo "... Two FIFO pileups prepared for $sample"
fi


# 04. VARIANT CALLING (SNP and INDELs)
# -------------------------------
echo "... Variant call for $sample"
java -Xmx8g -jar $varscn pileup2snp $fname --min-coverage $minCov --min-reads2 $minVar --min-avg-qual $minQ --p-value 1 -output-vcf 1 > "$sample.varscan.all.vcf"
rm $fname

echo $(pwd)

echo "... INDELs call for $sample"
java -Xmx8g -jar $varscn pileup2indel $fname_indel --min-coverage $minCov --min-reads2 $minVar --min-avg-qual $minQ --p-value 1 -output-vcf 1 > "$sample.varscan.indel.all.vcf"
rm $fname_indel

# 05. COMPUTE ON TARGET OF VCF
# ----------------------------
echo "... intersect with the target (SNP)"
head -n1 "$sample.varscan.all.vcf" > "$sample.varscan.all.BP.OT.vcf"
awk 'NR > 1 {print $1"\t"$2-1"\t"$2"\t"$0}' "$sample.varscan.all.vcf" | intersectBed -wa -a - -b $ftarget | cut -f 4-  | awk '{print $0}' >> "$sample.varscan.all.BP.OT.vcf"

echo "... intersect with the target (Indels)"
head -n1 "$sample.varscan.indel.all.vcf" > "$sample.varscan.indel.all.BP.OT.vcf"
awk 'NR > 1 {print $1"\t"$2-1"\t"$2"\t"$0}' "$sample.varscan.indel.all.vcf" | intersectBed -wa -a - -b $ftarget | cut -f 4-  | awk '{print $0}' >> "$sample.varscan.indel.all.BP.OT.vcf"

# 06. ANNOTATE WITH ANNOVAR
# --------------------------

# Annotate SNP
get_annovar $annovar_path $sample "$root" "$root/annovar_out"
mv "./annovar_out/$sample.hg19_multianno.txt" "$sample.annotation.BP.OT.txt"
rm -r ./annovar_out

# Annotate INDELs
get_annovar_indels $annovar_path $sample "$root" "$root/annovar_out"
mv "./annovar_out/$sample.hg19_multianno.txt" "$sample.annotation.indel.BP.OT.txt"
rm -r ./annovar_out



# 08. CALL DAMGING SNP
# ---------------------
get_dmg $sample "$root" "$root/$sample.dmg" $dbsnfs


















