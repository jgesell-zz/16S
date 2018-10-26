#!/usr/bin/env perl
#
use warnings;
use strict;



my $raw = $ARGV[0];
my $merged = $ARGV[1];
my $mapped = $ARGV[2];





my $reads = {};

my @raw;
open STATSMERGED, ">tempmerge";
open STATSMAPPED, ">temppick";
open STATSRAW, ">tempraw";

open IN, "$raw";
while (my $line = <IN>) {
	chomp $line;
	$line =~ s/^\s+//g;
	my @parts = split /\s+/, $line;
	$reads->{$parts[0]}->{'raw'} = $parts[1];
	print STATSRAW "$parts[1]\n";
}
close STATSRAW;
close IN;



open IN, "$merged";
while (my $line = <IN>) {
	chomp $line;
	$line =~ s/^\s+//g;
	my @parts = split /\s+/, $line;
#	print "Merged sample = $parts[0]\n";
	$reads->{$parts[1]}->{'merged'} = $parts[0];
	print STATSMERGED "$parts[0]\n";
}
foreach my $sample (keys %{$reads}) {
	unless (exists $reads->{$sample}->{'merged'}) {
		$reads->{$sample}->{'merged'} = 0
	}
	unless (exists $reads->{$sample}->{'raw'}) {
		$reads->{$sample}->{'raw'} = 0
	}


}
close STATSMERGED;
close IN;

open IN, "$mapped";
while (my $line = <IN>) {
        chomp $line;
	$line =~ s/^\s+//g;
        my @parts = split /\t/, $line;
#        print "Mapped sample =  $parts[1]\n";
	$reads->{$parts[0]}->{'mapped'} = $parts[1];
        print STATSMAPPED "$parts[1]\n";
}
foreach my $sample (keys %{$reads}) {
        unless (exists $reads->{$sample}->{'mapped'}) {
                $reads->{$sample}->{'mapped'} = 0
        }
        unless (exists $reads->{$sample}->{'raw'}) {
                $reads->{$sample}->{'raw'} = 0
        }


}
close STATSMAPPED;
close IN;

print "SampleID\tRaw\tMapped\tUnmapped\n";
foreach my $sample (sort {($reads->{$a}->{'mapped'} + ($reads->{$a}->{'raw'} - $reads->{$a}->{'mapped'})) <=> ($reads->{$b}->{'mapped'} + ($reads->{$b}->{'raw'} - $reads->{$b}->{'mapped'}))} keys %{$reads}) {
	my $mapped = $reads->{$sample}->{'mapped'};
	my $diff = $reads->{$sample}->{'merged'} - $mapped;
#	my $cmd = join(' cat ',$ENV{TMPDIR},'/',$sample,'.1.fq | wc -l');
#	my $readcount = `cat $ENV{TMPDIR}/$sample.1.fq | wc -l`/4;
	my $readcount = $reads ->{$sample}->{'raw'};
#	print STATSRAW "$readcount\n";
	print "$sample\t$readcount\t$mapped\t$diff\n";
}

print "\n\n";
print "Raw Stats\n";
my $capture = `cat tempraw | ~mcwong/listStats.pl`;
print "$capture";
`rm tempraw`;
print "\n\n";
print "Merged Stats:\n";
$capture = `cat tempmerge | ~mcwong/listStats.pl`;
print "$capture";
`rm tempmerge`;
print "\n\n";
print "Mapped Stats:\n";
$capture = `cat temppick | ~mcwong/listStats.pl`;
print "$capture";
`rm temppick`;
