#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use IO::File;

###loop through directory

#print "Analyzing...\n";
my $flag = 0;
my $v2;
my @content;

	my $fh = IO("fixedGC.out");
	my $w = O("fixedGC.out.gen");
	while(my $line = $fh->getline()) {
		chomp($line);		
		if($line=~/Sampled/) {
			chomp($line);
			$v2 = $line;
			$line =~s/^>\s+\w+\s+//;
			$line =~s/\s+.*$//;
			if($line == 0) {
				$flag = 0;
			} else {
				$v2 =~s/^>\s+\w+\s+\d+\s+.*(?=\d)//;
				$v2 =~s/\s.*$//;
				#need to add push
				push @content, "$v2\n";
				$flag = 1;
			}
		}	
		if($flag == 1) {
			if($line =~/\w+\s+/) {
				$line=~s/\s+//;
				push @content,"$line\n";
				$line=~s/\s+\-?\d+\.?\d+//;
				print $w uc $line."\n";
			}
		}
	}
	$w->close();
	$fh->close();
	my $wf = O("fixedGC.more");
		foreach(@content) {
			print $wf uc $_ ;
		}			
	$wf->close();

#print "Parsing Completed...\n";
#print "Please check ur current directory for ur files\n";
#print "tip to find ur files, enter this in terminal\n";
#print "ls * | grep *.in\n";

sub IO() {
	my $file = shift;
	return new IO::File($file,"r") or die "File not found: $!\n";
}
sub O() {
	my $file = shift;
	return new IO::File($file,"w+") or die "Cant write...: $!\n";
}
