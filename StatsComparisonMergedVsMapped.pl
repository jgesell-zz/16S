#!/usr/bin/env perl
#
use warnings;
use strict;



my $raw = $ARGV[0];
my $picked = $ARGV[1];





my $reads = {};
open IN, "$raw";

my @raw;
open STATSMERGED, ">tempmerge";
open STATSPICK, ">temppick";

while (my $line = <IN>) {
	chomp $line;
	$line =~ s/^\s+//g;
	my @parts = split /\s+/, $line;
	$reads->{$parts[1]}->{'raw'} = $parts[0];
	print STATSMERGED "$parts[0]\n";
}
close STATSMERGED;
close IN;



open IN, "$picked";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
#	print "sample = $parts[0]\n";
	$reads->{$parts[0]}->{'picked'} = $parts[1];
	print STATSPICK "$parts[1]\n";
}
foreach my $sample (keys %{$reads}) {
	unless (exists $reads->{$sample}->{'picked'}) {
		$reads->{$sample}->{'picked'} = 0
	}
	unless (exists $reads->{$sample}->{'raw'}) {
		$reads->{$sample}->{'raw'} = 0
	}


}


close STATSPICK;
print "SampleID\tRaw\tMapped\tUnmapped\n";
open STATSRAW, ">tempraw";
foreach my $sample (sort {($reads->{$a}->{'picked'} + ($reads->{$a}->{'raw'} - $reads->{$a}->{'picked'})) <=> ($reads->{$b}->{'picked'} + ($reads->{$b}->{'raw'} - $reads->{$b}->{'picked'}))} keys %{$reads}) {
	my $picked = $reads->{$sample}->{'picked'};
	my $diff = $reads->{$sample}->{'raw'} - $picked;
	my $cmd = join(' cat ',$ENV{TMPDIR},'/',$sample,'.1.fq | wc -l');
	my $readcount = `cat $ENV{TMPDIR}/$sample.1.fq | wc -l`/4;
	print STATSRAW "$readcount\n";
	print "$sample\t$readcount\t$picked\t$diff\n";
}
close STATSRAW;

print "\n\n";
print "Raw Stats\n";
my $capture = `cat tempraw | ~mcwong/listStats.pl`;
print "$capture";
`rm tempraw`;
print "\n\n";
print "Merged Stats:\n";
my $capture = `cat tempmerge | ~mcwong/listStats.pl`;
print "$capture";
`rm tempmerge`;
print "\n\n";
print "Mapped Stats:\n";
$capture = `cat temppick | ~mcwong/listStats.pl`;
print "$capture";
`rm temppick`;