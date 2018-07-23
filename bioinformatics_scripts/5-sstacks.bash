
samp=""
folder="ustacks-out-M3m3mls3"

while read name
do
  	samp+="-s ./${folder}/${name}_merged "
done < fastqlist.txt

sstacks -b 1 -c ./${folder}/batch_1 $samp -o ./${folder}/ -p 12

