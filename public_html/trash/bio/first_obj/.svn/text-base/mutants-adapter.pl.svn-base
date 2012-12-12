#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use IO::File;

###loop through directory

#print "Analyzing...\n";
my $content;
my $flag = 0;
$content = {};
	my $fh = IO("rmutants.out");
#	my $fh = IO("out.sample");
	my $w = O("rmutants.out.gen");
	while(my $line = $fh->getline()) {
		chomp($line);		
		if($line=~/superoptimal/) {
			$line = $fh->getline();
			chomp($line);
			$line=~s/\s+\(\-?\d+\.?\d+\)//;
			print $w uc $line."\n";
		}	
	}
	$w->close();
	$fh->close();

#print "Parsing Completed...\n";
#print "Please check ur current directory for ur files\n";
#print "tip to find ur files, enter this in terminal\n";
#print "ls * | grep *.in\n";

sub IO() {
	my $file = shift;
	return new IO::File($file,"r") or die "File not found: $!\n";
	my $w = O("rmutants.out.gen");
	while(my $line = $fh->getline()) {
		chomp($line);		
		if($line=~/superoptimal/) {
			$line = $fh->getline();
			chomp($line);
			$line=~s/\s+\(\-?\d+\.?\d+\)//;
			print $w uc $line."\n";
		}	
	}
#	$w->close();
	$fh->close();
}
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
