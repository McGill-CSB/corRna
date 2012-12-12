#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use IO::File;
use Getopt::Std;

our ($opt_o);
getopt('o');
###loop through directory
	my @files = <*.ps>;

#print "All files:: \n" . Dumper(@files) . "\n";
print "Analyzing...\n";
my $content;
$content = {};
my $fh;
foreach(@files) {
	if($opt_o) {
		if($_=~/^$opt_o/) {
		} else {
			next;
		}
	}
			
	 $fh = IO($_);
	my $w = O($_.".in");
	while(my $line = $fh->getline()) {
		chomp($line);		
		if($line=~/\/sequence/) {
			$line = $fh->getline();
			chomp($line);
			$line =~s/\\//;
			print $w $line ."\n";
		}	
		if($line=~/ubox/) {
			if($line=~/sqrt/ || $line =~/\/ubox/) {
			} else {
				print $w $line . "\n";
			}
		}
	}
	$w->close();
	$fh->close();
} 

print "Parsing Completed...\n";
print "Please check ur current directory for ur files\n";
print "tip to find ur files, enter this in terminal\n";
print "ls * | grep *.in\n";

sub IO() {
	my $file = shift;
	return new IO::File($file,"r") or die "File not found: $!\n";
}
sub O() {
	my $file = shift;
	return new IO::File($file,"w+") or die "Cant write...: $!\n";
}
