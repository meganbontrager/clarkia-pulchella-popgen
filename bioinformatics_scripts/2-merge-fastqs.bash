# piling up paired ends together since stacks doesn't have a means for dealing with them
# treat as independent loci for now and check for linkage later

while read name
do
cat clean-reads-trimmed/${name}*.fq.gz > merged-reads/${name}_merged.fq.gz
done < fastqlist.txt

