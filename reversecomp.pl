#!/usr/bin/env perl


use warnings;

use strict;

my $in = "";

if ($ARGV[0]) {
	$in = $ARGV[0];
} else {
	$in = "-";
}

open IN, "$in";


while (<IN>) {
	my $line = $_;
	chomp $line;
	$line = uc($line);
	$line =~ tr/ACGT UYRS WKMB DHVN/TGCA ARYS WMKV HDBN/;
	$line = reverse($line);
	print "$line\n";
}


