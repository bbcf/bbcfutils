#!/bin/sh
arguments=$*
echo "arguments: "${arguments}

#mergeBigWig.sh -i file1,file2,file3...filen -o outputfile -a assembly
#./mergeBigWig.sh -i density_file_CHAN1_2nM_A_filtered_sorted.bw,density_file_CHAN1_2nM_B_filtered_sorted.bw -o fileres -a hg19

while getopts i:o:a: name   ## les options acceptant un paramètres sont suivies de ":"
  do
    case $name in
        i)
            bwFiles_flag=1
            bwFiles="$OPTARG"
            ;;
        o)
            outputFile_flag=1
            outputFile="$OPTARG"
            ;;
	a)
	    assembly_flag=1
	    assembly="$OPTARG"
	    ;;
        ?)
            printf "Usage: %s: [-i bigWig files (comma separated)] [-o output filename] [-a assembly (if provided result will be in a bigWig format)] args\n" $0
            exit 2
            ;;
    esac
done


infiles=""
for x in `echo $bwFiles|awk 'BEGIN{FS=","}{for(i=1;i<=NF;i++){print $i}}'`; do
#echo "Will transform file $x to bedGraph file"
`bigWigToWig $x $x".tmp"`
`grep -v '#bedGraph' $x".tmp" > $x".bedGraph"`
infiles=$infiles$x".bedGraph,"
rm $x".tmp"
done

infiles=`echo $infiles | sed 's/,/ /g'`
tmpfile=$outputFile".tmp"
`cat $infiles  | awk '{for(i=$2;i<$3;i++){print $1"\t"i"\t"$4}}' | sort -k 1,1 -k2,2n | groupBy -i stdin -grp 1,2 -c 3 -o mean > $tmpfile`
awk '{if(NR==1){ch=$1;s=$2;e=$2;curScore=$3} if($1==ch && $3==curScore){e=$2}else{e++;print ch"\t"s"\t"e"\t"curScore; ch=$1;s=$2;e=$2;curScore=$3}} END{e++;print ch"\t"s"\t"e"\t"curScore}' $tmpfile  > $outputFile 
if [[ ! -z "${assembly_flag}" ]]
then
	echo "will convert bedGraph to bw (assembly="$assembly")"
	`genrep4humans.py -l -a $assembly|cut -f2,3 > chromList.txt` 
	`wigToBigWig $outputFile chromList.txt $tmpfile`
	mv $tmpfile $outputFile	
fi

rm $tmpfile
rm $infiles
