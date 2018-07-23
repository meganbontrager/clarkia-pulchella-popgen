# running rxstacks and then redoing cstacks and sstacks

rxstacks -b 1 -P ustacks-out-M3m3mls3/ -o corr-out-M3m3mls3/ --lnl_filter --lnl_lim -10.0 --conf_filter --conf_lim 0.75 --prune_haplo

samp=""
folder="corr-out-M3m3mls3"

while read name
do
    samp+="-s ${folder}/${name}_merged "
done < fastqlist.txt

cstacks -b 1 -o ./${folder} -n 3 -p 12 $samp

samp=""
folder="corr-out-M3m3mls3"

while read name
do
  	samp+="-s ./${folder}/${name}_merged "
done < fastqlist.txt

sstacks -b 1 -c ./${folder}/batch_1 $samp -o ./${folder}/ -p 12
