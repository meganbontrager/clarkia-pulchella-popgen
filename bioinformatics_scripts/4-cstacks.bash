
samp=""
folder="ustacks-out-M3m3mls3"

# make a long concatenated list of samples to create the catalog
# later read in a catchen publication that it's better to use fewer samples to build the catalog
# so if redoing, amend that

while read name
do 
    samp+="-s ${folder}/${name}_merged "
done < fastqlist.txt

cstacks -b 1 -o ./${folder} -n 3 -p 12 $samp

