#!/usr/bin/perl -w
use strict;

die "perl $0 <in><cut-off>\n" unless(@ARGV==2);

open I, "$ARGV[0]" or die;
my %raw;
my ($abund_s,$uniq_s,$num) = (0,0,0);
while(<I>)
{
	chomp;
	my $id = $_;
	chomp(my $seq=<I>);
	my ($abund,$uniq) = (split /:/,$id)[3,4];
	$raw{$id} = [($abund,$uniq,$seq)];
	$abund_s += $abund;
	$uniq_s += $uniq;
	$num++;
}
close I;

my $cutoff=$ARGV[1];
for(keys %raw){
	my ($abund,$uniq,$seq) = @{$raw{$_}};
	my $abund_r = $abund/$abund_s*100*$num;
	my $uniq_r = $uniq/$uniq_s*100*$num;
	if($abund_r>$cutoff && $uniq_r>$cutoff){
		print "$_\n$seq\n";
	}
	
}
