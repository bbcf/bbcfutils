#!/bin/sh
arguments=$*
echo "arguments: "${arguments}

#mergeBigWig.sh -i file1,file2,file3...filen -o outputfile
#./mergeBigWig.sh -i density_file_CHAN1_2nM_A_filtered_sorted.bw,density_file_CHAN1_2nM_B_filtered_sorted.bw -o fileres

while getopts i:o: name   ## les options acceptant un paramètres sont suivies de ":"
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
        ?)
            printf "Usage: %s: [-i bigWig files (comma separated)] [-o output filename] args\n" $0
            exit 2
            ;;
    esac
done

infiles=""
for x in `echo $bwFiles|awk 'BEGIN{FS=","}{for(i=1;i<=NF;i++){print $i}}'`; do
#echo "Will transform file $x to bedGraph file"
`bigWigToWig $x $x".tmp"`
`grep -v '#bedGraph' $x".tmp" > $x".bedGraph"`
#`awk '{if($0 !~ /^#bedGraph/){print > FILENAME"_"$1".bedGraph"}}' $x".tmp"`
infiles=$infiles$x".bedGraph,"
rm $x".tmp"
done

#allChr=for x in `cut -f1 density_file_CHAN1_2nM_B_filtered_sorted.bw.bedGraph|uniq`;do echo "curChr="$x;done


#./mergeBigWig.sh -i density_file_CHAN1_2nM_A_filtered_sorted.bw,density_file_CHAN1_2nM_B_filtered_sorted.bw -o fileres

infiles=`echo $infiles | sed 's/,/ /g'`
#echo "infiles="$infiles
tmpfile=$outputFile".tmp"
#echo "tmpfile="$tmpfile
`cat $infiles  | awk '{for(i=$2;i<=$3;i++){print $1"\t"i"\t"$4}}' | sort -k 1,1 -k2,2n | groupBy -i stdin -grp 1,2 -c 3 -o mean > $tmpfile`
awk '{if(NR==1){ch=$1;s=$2;e=$2;curScore=$3} if($1==ch && $3==curScore){e=$2;new=0}else{print ch"\t"s"\t"e"\t"curScore; ch=$1;s=$2;e=$2;curScore=$3;new=1}} END{print ch"\t"s"\t"e"\t"curScore}' $tmpfile  > $outputFile 
rm $tmpfile
rm $infiles
