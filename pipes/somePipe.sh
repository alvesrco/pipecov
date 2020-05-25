SAMPLE_R1=$1
SAMPLE_R2=$2
RESULTS=$3
REFERENCE=$4

#*************************************************

IFS='/' tokens=( $SAMPLE_R1 )
SAMPLE_NAME=${tokens[1]}
DIR=$RESULTS/$SAMPLE_NAME
IFS=' '
[ -d $DIR ] || mkdir $DIR

echo "Current sample: $SAMPLE_NAME"

#limpiar fastq con fastp o prinseq
echo "Cleaning FASTQ...."
IN1=$SAMPLE_R1
IN2=$SAMPLE_R2
OUT1=$DIR/$SAMPLE_NAME\_1_fastp.fastq
OUT2=$DIR/$SAMPLE_NAME\_2_fastp.fastq
LOG=$DIR/$SAMPLE_NAME\_fastp.log
fastp --thread 8 -i $IN1 -I $IN2 -o $OUT1 -O $OUT2 > $LOG 2>&1
fastqc $OUT1 $OUT2 -o $DIR -q

#remover duplicados con bbmap
echo "Removing duplicates..."
IN1=$OUT1
IN2=$OUT2
OUT1=$DIR/$SAMPLE_NAME\_1_dedupped.fastq
OUT2=$DIR/$SAMPLE_NAME\_2_dedupped.fastq
LOG1=$DIR/$SAMPLE_NAME\_1_dedupped.log
LOG2=$DIR/$SAMPLE_NAME\_2_dedupped.log
dedupe.sh in=$IN1 out=$OUT1 > $LOG1 2>&1
dedupe.sh in=$IN2 out=$OUT2 > $LOG2 2>&1
fastqc $OUT1 $OUT2 -o $DIR -q

#deconseq para sacar contaminación humana dejando lo que esté en sarscov2
echo "Decontaminating..."
IN1=$OUT1
IN2=$OUT2
OUT1=$DIR/$SAMPLE_NAME\_1_clean.fq
OUT2=$DIR/$SAMPLE_NAME\_2_clean.fq
deconseq.pl -f $IN1 -dbs hsref -dbs_retain sarscov2 -id $SAMPLE_NAME\_1 -out_dir $DIR
deconseq.pl -f $IN2 -dbs hsref -dbs_retain sarscov2 -id $SAMPLE_NAME\_2 -out_dir $DIR
fastqc $OUT1 $OUT2 -o $DIR -q

#reparar los fastq para que coincidan la cantidad y orden de los reads
echo "Re-pairing FASTQ..."
IN1=$OUT1
IN2=$OUT2
OUT1=$DIR/$SAMPLE_NAME\_1_repaired.fastq
OUT2=$DIR/$SAMPLE_NAME\_2_repaired.fastq
SINGLE=$DIR/$SAMPLE_NAME\_repaired_single.fastq
repair.sh in1=$IN1 in2=$IN2 out1=$OUT1 out2=$OUT2 outsingle=$SINGLE

#Mapeo
echo "Mapping..."
IN1=$OUT1
IN2=$OUT2
OUT1=$DIR/$SAMPLE_NAME.sam
OUT2=$DIR/$SAMPLE_NAME.bam
LOG=$DIR/$SAMPLE_NAME\_bwa.log
STAT=$DIR/$SAMPLE_NAME\_flagstat.txt
bwa mem -t 8 $REFERENCE $IN1 $IN2 > $OUT1 2>$LOG
samtools sort $OUT1 -o $OUT2
samtools index $OUT2
samtools flagstat $OUT2 > $STAT

#Llamado de variantes
echo "Calling variants..."
IN1=$OUT2
OUT1=$DIR/$SAMPLE_NAME\_mpileup.vcf.gz
OUT2=$DIR/$SAMPLE_NAME\_call.vcf.gz
bcftools mpileup -f $REFERENCE $IN1  -Oz > $OUT1
bcftools call -mv -Oz $OUT1 -o $OUT2

#Obtención de secuencia consenso
echo "Getting consensus sequence..."
IN1=$OUT2
OUT1=$DIR/$SAMPLE_NAME\_consensus.fasta
bcftools index $IN1
bcftools consensus -f $REFERENCE -H A -o $OUT1 $IN1

#Assembly
echo "Getting assembly..."
IN1=$DIR/$SAMPLE_NAME\_1_repaired.fastq
IN2=$DIR/$SAMPLE_NAME\_2_repaired.fastq
megahit -1 $IN1 -2 $IN2 -o assembly