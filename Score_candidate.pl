#!/usr/bin/perl -w
use strict;

die "perl $0 <in><raw.log><combined.log><combined.file>\n" unless(@ARGV==4);
#
# print out format:
# print "#id title:\n#Gene-infer-flag:Score:total_seq_support:uniq_seq_support:combined_seq_num\n";

my ($in, $rawlog, $combined_log, $combined_f) = @ARGV;

my %score;

my ($max_trim5_abund, $max_trim5_uniq, $max_more5, $max_dep_abund,$max_dep_uniq , $max_combined_num) = (0,0,0,0,0,0);
#----------	1. Get the score of depth rate (abundance && unique) 
open C, "$combined_f" or die;

my %raw;
my ($abund_s,$uniq_s,$num) = (0,0,0);

while(<C>)
{
	chomp;
	s/>//;
	my @t = split /:/,$_;
        my $id = $_;
        chomp(my $seq=<C>);
        my ($abund,$uniq) = (split /:/,$id)[3,4];
	my $id_trim = "$t[0]:$t[1]";;
        $raw{$id_trim} = [($abund,$uniq,$seq)];
        $abund_s += $abund;
        $uniq_s += $uniq;
        $num++;	
}
close C;

for(keys %raw){
        my ($abund,$uniq,$seq) = @{$raw{$_}};
        my $abund_r = $abund/$abund_s*100*$num;
        my $uniq_r = $uniq/$uniq_s*100*$num;
	$max_dep_abund = $abund_r if($max_dep_abund<$abund_r);
	$max_dep_uniq = $uniq_r if($max_dep_uniq<$uniq_r);
	$score{$_}->[0] = $abund_r;
	$score{$_}->[1] = $uniq_r;
}

#-----------------	2. Get the score of trimed ratios(abundance && unique) and more5 ratio

my %raw_back;
open R, "$rawlog" or die;
while(<R>)
{
	chomp;
	s/^>//;
	my @line = split;
	my @t = split /:/,$line[0];
	my $id = "$t[0]:$t[1]";
	$raw_back{$id} = "$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[-1]";
}
close R;

my %combine_back;
open L, "$combined_log" or die;
while(<L>)
{
	chomp;
	my @line = split;
	if($#line==1){
		push @{$combine_back{$line[0]}}, $line[1];
	}else{
		$combine_back{$line[0]}->[0] = $line[0];
	}
}
close L;

for my $id1(keys %combine_back)
{
	my ($abund1, $uniq1, $abund2,$uniq2 , $more5) = split /\t/,$raw_back{$id1};
	my $num = scalar(@{$combine_back{$id1}}) +1;
	if($combine_back{$id1}->[0] eq $id1){
		$num = 1;
	}else{
		for my $id2(@{$combine_back{$id1}})
		{
		my @line = split /\t/,$raw_back{$id2};
		$abund1 += $line[0];
		$uniq1 += $line[1];
		$abund2 += $line[2];
		$uniq2 += $line[3];
		$more5 += $line[4];
		}
	}
	my $abund_s = $abund2/$abund1;
	my $uniq_s = $uniq2/$uniq1;
	my $more5_s = $more5/$num;
	
	$max_trim5_abund = $abund_s if($max_trim5_abund<$abund_s);
	$max_trim5_uniq = $uniq_s if($max_trim5_uniq<$uniq_s);
	$max_more5 = $more5_s if($max_more5<$more5_s);

	$score{$id1}->[2] = $abund_s ;
	$score{$id1}->[3] = $uniq_s ;
	$score{$id1}->[4] = $more5_s;
}

# ---------------	3 Get Combined_seq && Major/Minor scores

open I, "$in" or die;
while(<I>)
{
	s/>//;
	my $id = $_;
	chomp(my $seq = <I>);
	my @t = split /:/,$id;
#	my $id_trim = "$t[0]:$t[1]";
	$max_combined_num = $t[2] if($max_combined_num < $t[2]);
}
close I;

open I, "$in" or die;
my $flag = 0;
#print "#id title:\n#Gene-infer-flag Score total_seq_support uniq_seq_support combined_seq_num\n";
while(<I>)
{
	$flag++;
	s/>//;
	my $id = $_;
	chomp(my $seq = <I>);
	my @t = split /:/,$id;
	my $gene = (split /_/,$t[0])[0];
	my $id_trim = "$t[0]:$t[1]";
	my $major_minor_s = 15;
	$major_minor_s = 10 if($t[-1] == 2);
	my $S = sqrt($score{$id_trim}->[2]/$max_trim5_abund)*20+sqrt($score{$id_trim}->[3]/$max_trim5_uniq)*20+sqrt($score{$id_trim}->[4]/$max_more5)*15+sqrt($score{$id_trim}->[0]/$max_dep_abund)*10+sqrt($score{$id_trim}->[1]/$max_dep_uniq)*10+sqrt($t[2]/$max_combined_num)*15+sqrt($major_minor_s);
	$S = int $S;
	print ">$gene-infer-$flag:$S:$t[3]:$t[4]:$t[2]\n$seq\n";
}
close I;

