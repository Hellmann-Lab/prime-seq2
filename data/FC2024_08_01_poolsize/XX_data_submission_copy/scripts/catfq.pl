#!/usr/bin/perl
# LMU Munich. AG Enard
# A script to filter reads based on Barcode base quality.
# Author: Swati Parekh
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de
# Modified by Daniel Richter to include i5 reads of new dual-indexed layout!

if(@ARGV != 5)
{
print 
"\n#####################################################################################
Usage: perl $0 <i7-Read.fq.gz> <barcode-Read.fq.gz> <output.fq> <threads>\n
Explanation of parameter:

i7-Read.fq.gz	- Input i7 reads fastq file name.
i5-Read.fq.gz	- Input i5 reads fastq file name.
barcode-Read.fq.gz- Input trimmed R1 reads fastq file name (needs to contain only sample barcode)
output.fq	- Output file name. pigz will put the .gz only provide the base name.
threads		- number of processors to zip.
Please drop your suggestions and clarifications to <parekh\@bio.lmu.de>\n
######################################################################################\n\n";
exit;
}

$i7read=$ARGV[0];
$i5read=$ARGV[1];
$bcread=$ARGV[2];
$bcreadoutfull = $ARGV[3];
$threads=$ARGV[4];

if ($bcread =~ /\.gz$/) {
open BCF, '-|', 'gzip', '-dc', $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open I7F, '-|', 'gzip', '-dc', $i7read || die "Couldn't open file $i7read. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open I5F, '-|', 'gzip', '-dc', $i5read || die "Couldn't open file $i5read. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}
else {
open BCF, "<", $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open I7F, "<", $i7read || die "Couldn't open file $i7read. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open I5F, "<", $i5read || die "Couldn't open file $i5read. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}

open BCOUTFULL, ">", $bcreadoutfull || die "Couldn't open file $bcreadoutfull to write\n\n";;

#these counters are actually not used within this script; commented out for now...
#$count=0;
$total=0;
#$filtered=0;

while(<BCF>){
	$total++;
	$brid=$_;
	$brseq=<BCF>;
	chomp($brseq);
	$bqid=<BCF>;
	$bqseq=<BCF>;
	chomp($bqseq);
	
	$i7rid=<I7F>;
	$i7rseq=<I7F>;
	chomp($i7rseq);
	$i7qid=<I7F>;
	$i7qseq=<I7F>;
	chomp($i7qseq);
	
	$i5rid=<I5F>;
	$i5rseq=<I5F>;
	chomp($i5rseq);
	$i5qid=<I5F>;
	$i5qseq=<I5F>;
	chomp($i5qseq);
	
	$seq=$i7rseq.$i5rseq.$brseq;
	$qseq=$i7qseq.$i5qseq.$bqseq;
	print BCOUTFULL $brid,$seq,"\n",$bqid,$qseq,"\n";
}
close BCF;
close I7F;
close I5F;
close BCOUTFULL;
`pigz -f -p $threads $bcreadoutfull`;

print "Finished! Merged number of reads: ".$total;

