#!/usr/bin/perl -w
use Parallel::ForkManager;
use Data::Dumper;

# parameter file name as commandline argument
my $paramfile = $ARGV[0]; 

# read in param file, PARAM is the file handle, and create a hash of all the parameters.
# split by => and make the first column, the key and the second column, the value. 
open PARAM, $paramfile or die print $!;
my %param;
while(<PARAM>)
{
	chomp;
	@r = split('=>');
	#print "@r";
	$param{$r[0]}=$r[1];
}

$param{'ALNDIR'} = $param{'PROJECTNAME'}.'/GSNAP';

# make directory with projectname and subdirectories gsnap, cufflinks and genecounts
# system("mkdir $param{'PROJECTNAME'}") unless (-d $param{'PROJECTNAME'});
# system("mkdir $param{'ALNDIR'}") unless (-d $param{'ALNDIR'});

# open file with filename which is value of the key FASTALIST
open FILE, $param{'FASTALIST'} or die;

# open file for writing
my $filename = $param{'ALNDIR'}."/report.csv";
open my $fh, ">", $filename or die("Could not open file. $!");

# write header to report
@header = ("Sample",'Library','ReadsMapped','ReadsUnmapped','TotalReads','PercReadsMapped','PercReadsUnmapped');
$header = join(',',@header);
$header = $header."\n";
print $fh $header;

# create an array samples
my @samples;

# splitting each line based on ',' and storing in an array @r
# pushing the reference of this array in another array @samples
while(<FILE>){
	chomp;
	my @r = split(',');
	push(@samples,\@r);
}

my $start_run = time();

# run parallel jobs
my $pm=new Parallel::ForkManager($param{'JOBS'});

foreach (@samples)
{	
	$pm->start and next;		
        
        	# use cutadapt
        	my $cutadapt = "/fujfs/d3/home/komalr/tools/cutadapt-master/bin/cutadapt --untrimmed-output=$param{'ALNDIR'}/".$_->[1]."_untrimmed.fa -a short_adpt=$param{'short_adapter'} -g long_adpt=$param{'long_adapter'} --times=2 --error-rate=$param{'ERROR_RATE'} --output=$param{'ALNDIR'}/".$_->[1]."_{name}_trimmed.fa --info-file=$param{'ALNDIR'}/".$_->[1]."_info.txt $param{'FASTQDIR'}/$_->[0]";
        	print $cutadapt,"\n";
            	system($cutadapt);
        
        	# keep first 20 bp from long-clipped file, automatically clip short adapters that are not clipped by cutadapt because of error rate
        	my $trimends = "sed -e '2~2s/^\\(.\\{20\\}\\).*/\\1/' $param{'ALNDIR'}/".$_->[1]."_long_adpt_trimmed.fa > $param{'ALNDIR'}/".$_->[1]."_long_adpt_trimmed-1.fa";
        	print $trimends,"\n";
            	system($trimends);
        
        	# remove tmp file
        	my $movefile = "mv $param{'ALNDIR'}/".$_->[1]."_long_adpt_trimmed-1.fa $param{'ALNDIR'}/".$_->[1]."_long_adpt_trimmed.fa";
        	print $movefile,"\n";
            	system($movefile);
        
        	# concatenate files
        	my $concat = "cat $param{'ALNDIR'}/".$_->[1]."*trimmed.fa | gzip - > $param{'ALNDIR'}/".$_->[1]."_processed.fa.gz";
        	print $concat,"\n";
            	system($concat);
        
        	# run gsnap
        	my $gsnapcmd = "gsnap -d $_->[2] --gunzip --batch=5 -m 2 --nthreads=$param{'THREADS'} --format=sam $param{'ALNDIR'}/".$_->[1]."_processed.fa.gz 2> $param{'ALNDIR'}/$_->[1].out | samtools view -bS - > $param{'ALNDIR'}/$_->[1].bam";
        	print $gsnapcmd,"\n";
            	system($gsnapcmd);
            
            	# sort bam
            	my $samsort = "samtools sort -n $param{'ALNDIR'}/$_->[1].bam $param{'ALNDIR'}/".$_->[1]."_sorted";
            	print $samsort,"\n";
            	system($samsort);

            	# get mappings
            	my $mappings = "samtools view $param{'ALNDIR'}/".$_->[1]."_sorted.bam | grep MGLib | cut -f3 | sort | uniq -c | awk '{ t=\$1; \$1=\$2; \$2=t; print }' | sed 's/ /,/g' > $param{'ALNDIR'}/".$_->[1]."_UID_mappings.txt";
            	print $mappings,"\n";
            	system($mappings);
            
            	# get total reads
            	my $totalreads = `cat $param{'FASTQDIR'}/$_->[0] | wc -l`;
            	$totalreads = int($totalreads/2);
            	print $totalreads,"\n";
            
            	# get mapped reads
            	my $mapped = `samtools view $param{'ALNDIR'}/$_->[1]_sorted.bam | grep MG | cut -f1 | uniq | wc -l`;
            	$mapped = int($mapped);
        	print $mapped,"\n";
            	my $mappedperc = sprintf("%.2f", $mapped/$totalreads*100);
            
            	# get unmapped reads
            	my $unmapped = `samtools view $param{'ALNDIR'}/$_->[1]_sorted.bam | grep -v MG | cut -f1 | uniq | wc -l`;
            	$unmapped = int($unmapped);
            	print $unmapped,"\n";
            	my $unmappedperc = sprintf("%.2f", $unmapped/$totalreads*100);
            
            	# make list
            	my @list = ($_->[1],$_->[3],$mapped,$unmapped,$totalreads,$mappedperc,$unmappedperc);
            	my $list = join(',',@list);
            	$list = $list."\n";
            	print $fh $list;
            
		print "$_->[1] completed\n"; 
	
	$pm->finish;
}

$pm->wait_all_children;

close $fh; 

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";
