
# using stacks version 1.40
# beginning with raw reads, output of this script is the de-multiplexed reads that are available on SRA

process_radtags -i gzfastq -p ./raw-reads/lane1raw -P -o ./clean-reads-trimmed/ -b ./barcodes/barcodes-lane1.txt -t 93 -r -c -q -D --inline_inline --renz_1 pstI --renz_2 mspI -E phred33 --barcode_dist_1 1
mv clean-reads-trimmed/process_radtags.log clean-reads-trimmed/process_radtags_l1.log
process_radtags -i gzfastq -p ./raw-reads/lane2raw -P -o ./clean-reads-trimmed/ -b ./barcodes/barcodes-lane2.txt -t 93 -r -c -q -D --inline_inline --renz_1 pstI - renz_2 mspI -E phred33 --barcode_dist_1 1
mv clean-reads-trimmed/process_radtags.log clean-reads-trimmed/process_radtags_l2.log
