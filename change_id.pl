#!/usr/bin/perl -w
use strict;

die "perl $0 <in><filter.txt><out>\n" unless(@ARGV==3);

my %id_all;
open I, "$ARGV[1]" or die;
while(<I>)
{
	chomp;
	my @line = split;
	my @id = split /:/,$line[0];
	my $new = "$id[0]:$id[1]:$line[3]:$line[4]";
	$id_all{$line[0]} = $new;
}
close I;

open I, "$ARGV[0]" or die;
open O, ">$ARGV[2]" or die;
while(<I>)
{
	chomp;
	print O "$id_all{$_}\n";
	my $seq = <I>;
	print O "$seq";
}
close I;

