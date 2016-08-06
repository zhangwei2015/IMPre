#!/usr/bin/perl -w 
use strict;

die "perl $0 <in.list><out-dir><out-number>\n" unless(@ARGV==3);

my ($list , $out , $max) = @ARGV;
my %all;
my %id_all;
open I, "$list" or die;
while(<I>)
{
	chomp;
	next unless( -s $_);
	open N, "gzip -dc $_|" or die;
	chomp(my $id = <N>);
	my ($abund, $uniq);
	($id, $abund, $uniq) = split /:/,$id;
	if(exists $id_all{$id}){
		$id_all{$id}->[0] += $abund;
		$id_all{$id}->[1] += $uniq;
	}
	else{
		$id_all{$id}->[0] = $abund;
		$id_all{$id}->[1] = $uniq;
	}
	
	while(<N>){
		my $r = $_;
		my $s = <N>;
		if(exists $all{$id}){
			$all{$id} .= "$r$s";
		}else{
			$all{$id} = "$r$s";
		}
	}
	close N;
}
close I;

my $flag = 1;
for(sort {$id_all{$b}->[0] <=> $id_all{$a}->[0]} keys %id_all)
{
	last if($flag == $max);
	open O, "|gzip >$out/seed.$flag.seq.gz" or die;
	print O "$_:$id_all{$_}->[0]:$id_all{$_}->[1]\n";
	print O "$all{$_}";
	close O;
	$flag++;
}

