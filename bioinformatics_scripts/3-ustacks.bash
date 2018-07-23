# tried many different options for -m -M and --max_locus_stacks but this is what I went with in the end

mkdir ustacks-out-M3m3mls3
x=1
while read name
do
ustacks  -t gzfastq -f merged-reads/${name}_merged.fq.gz -o ustacks-out-M3m3mls3/ -i $x -m 3 -M 3 -H -r -d --max_locus_stacks 3 -p 12
x=$(($x+1))
done < fastqlist.txt
