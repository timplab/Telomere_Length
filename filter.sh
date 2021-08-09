#!/bin/bash 

while :
do
  case "$1" in 
    -i | --in ) # input Tidehunter
      in=$2
      shift 2
      ;;
    -r | --reads ) # fastq
      reads=$2
      shift 2
      ;;
    -s | --sequence ) #telomere sequence (optional : defaults to GTGTG)
      sequence=$2
      shift 2
      ;;
    -p | --seqkit_path ) #telomere sequence (optional : defaults to GTGTG)
      seqkit_path=$2
      shift 2
      ;;
    -h | --help ) # help message
      helpmsg=1
      shift 1
      ;;
    *) break
      ;;
  esac
done

if [ -z $in ];then
  errormsg="Error : Tidehunter output must be given"
  helpmsg=1
fi
if [ -z $reads ];then
  errormsg="Error : fastq must be given"
  helpmsg=1
fi
if [ ! -z $helpmsg ];then
  echo "./filter.sh [-i /Tidehunter/directory] [-s telomere_sequence [string]] /path/to/basecalled/sample
  Filter Tidehunter fastq for reads with telomere and adapter seqence
  ;-i/--in;;Tidehunter input file
  ;-r/--reads;;Fastq
  ;-s/--sequence;;Telomere sequence (default GTGTGTGGGTGTG)
  ;-p/--seqkit_path;;path to seqkit executable
  ;-h/--help;;show help message and exit"|\
    tr ";" "\t"
  echo $errormsg
  exit
fi
if [ -z $sequence ];then
  sequence=GTGTGTGGGTGTG
fi

mkdir -p ./tmp
out=./tmp

grep $sequence $in | cut -f1 | sed 's/^/@/' > ${out}/ID.tmp
#num=$(grep $sequence $in | cut -f1 | wc -l)
#echo $num reads with $sequence

##Generate a filtered .fastq with only reads with telomere
awk 'NR==FNR{a[$0];next}$1 in a {x=NR+3}(NR<=x){print} ' ${out}/ID.tmp ${reads} > ${out}/filtered.fastq
##Filter out reads < 5kb
#awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 5000) {print header, seq, qheader, qseq}}' < ${out}/telomere.fastq > ${out}/filtered.fastq
##Filter for reads with more than 10 As in the last 100 bp
$seqkit_path grep -s -R -100:-1 -r -p AAAAAAAAAA ${out}/filtered.fastq > ${out}/As.fastq
##Filter for reads with more then 10 Ts in the first 100 bp
$seqkit_path grep -s -R 1:100 -r -p TTTTTTTTTT ${out}/filtered.fastq > ${out}/Ts.fastq
##Combine the reads that are tailed
cat ${out}/As.fastq ${out}/Ts.fastq
rm -r $out
