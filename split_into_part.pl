#!/usr/bin/perl -w
use strict;

die "perl $0 <in.gz><num><out-dir><name>\n" unless(@ARGV==4);

my ($in,$num,$out,$name) = @ARGV;

my %all;
open I, "gzip -dc $in|" or die;
while(<I>)
{
	my $id = $_;
	chomp(my $seq= <I>);
	next if(length($seq) <=40);
	my $i = int rand($num);
	$all{"$id$seq"} = $i;
}
close I;


for(my $i=0 ; $i<$num; $i++)
{
	my %part;
	for (keys %all){
		$part{$_} = 1 if($all{$_}==$i);
	}
	&output($i,%part);
}


sub output
{
	my ($i,%h) = @_;
	$i++;
	open O, "|gzip >$out/$name.$i.gz" or die;
	for(keys %h)
	{
		print O "$_\n";
	}
}
