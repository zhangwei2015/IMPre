#!/usr/bin/perl -w
use strict;

die "perl $0 <*.gz_dir><*.fa> <miss_num><out><interval_filter><ave_r_filter><filter_more5-1><filter_more5-2><gene>\n" unless(@ARGV==9);

my ($dir,$fa,$interval_num ,$out , $interval_filter ,$ave_r_filter,$filter_more5_1,$filter_more5_2,$gene) = @ARGV;

#my $interval_filter = 1.5;# #(subseq of interval)/#(with $num_check) >1.5, both abundance and unique
#my $mini_abund = 50;# the subsequence with interval_cut
#my $ave_r_filter = 1; # the abundance_num/average_abund and unique_num/average_uniq >1%
# print out format:
# >V_100:1:7:7
# V_seed_name:flag:abundance:unique

open O, ">$out" or die;
#my $min_len = 30;

# --- read the germline sequence && in store
open I, "$ARGV[1]" or die;

my %germline;
my %backup;
while(<I>)
{
	chomp;
	my $raw_id = $_;
	my $flag = (split /[_:]/,$_)[1];
	chomp(my $seq = <I>);
	$backup{$raw_id} = $seq;

	my $subseq;
	if($gene eq "V"){
		$subseq= substr($seq,0,length($seq)-$interval_num);
	}else{
		$subseq= substr($seq,$interval_num);
	}
#	my $min_seed = substr($subseq,-$min_len);
	$germline{$flag}{$raw_id} = [( $subseq , $seq)];
}
close I;

#----	calulate the number of sequence surpporting the germline
#    supportint: part of sequence same with the germline is OK. for a V example as follow: two seq are OK
#    germline:	       ATGGACAGCCTGAGAGCCGAGGACACGGCTGTATATTAC
#	seq1:  actgacttATGGACAGCCTGAGAGCCGAGGACACGGCTGTATATTACgggg
#	seq2:                     GAGAGCCGAGGACACGGCTGTATATTACacttttttttttttttt


my %germline_stat;

my ($sum_sub_a,$sum_sub_u,$sum_seq_a,$sum_seq_u,$sum_more5) = (0,0,0,0,0);

for my $file (keys %germline)
{
	for my $id (keys %{$germline{$file}})
	{
		open I, "gzip -dc $dir/seed.$file.seq.gz|" or die;
		chomp(my $seed = <I>);
		$seed = (split /:/, $seed)[0];
		my $germ_seed_s = index($germline{$file}{$id}->[1],$seed) + 1;# get the position of seed in the germline
		
		@{$germline_stat{$id}}[0,1,2,3,4] = (0,0,0,0,0);
		my %more5_uniq;
		while(<I>)
		{
			chomp;
			my $id_seq = $_;
			chomp(my $seq = <I>);
			my ($start , $abund) = (split /[:-]/,$id_seq)[2,1];
			$germline_stat{$id}->[4] += $abund;
			$germline_stat{$id}->[5]++;
			my ($sub_germ,$germ) = ($germline{$file}{$id}->[0],$germline{$file}{$id}->[1]);
			if($gene eq "V" && $start < $germ_seed_s)# for V: only part of germline sequence in the seq, then cut the extra ones
			{
				$sub_germ = substr($germline{$file}{$id}->[0] , $germ_seed_s-$start);
				$germ = substr($germline{$file}{$id}->[1] , $germ_seed_s-$start);
			}

			if($seq =~ /$sub_germ/){
				$germline_stat{$id}->[2] += $abund;
				$germline_stat{$id}->[3]++;
				if($seq =~ /$germ/){
					$germline_stat{$id}->[0] += $abund;
					$germline_stat{$id}->[1]++;
					my $more5;
					$more5 = $gene eq "V" ?substr($',0,5):substr($`,-5,5);
					$more5_uniq{$more5} = 1;
				}
			}
		}	
		
		close I;
		$sum_seq_a += $germline_stat{$id}->[0];
		$sum_seq_u += $germline_stat{$id}->[1];
		$sum_sub_a += $germline_stat{$id}->[2];
		$sum_sub_u += $germline_stat{$id}->[3];
		my $more5_num = scalar(keys %more5_uniq);
		$germline_stat{$id}->[6] = $more5_num;
		$sum_more5 += $more5_num;
	}
}


# calculate rate && filter && output
# filter:
#	1. abund rate, uniq rate
#	2. average abund rate, average uniq rate

for my $id( keys %germline_stat)
{
	my ($abund_r ,$uniq_r) = (0,0);#calculate rate
	if($germline_stat{$id}->[0]){
		($abund_r ,$uniq_r) = ($germline_stat{$id}->[2]/$germline_stat{$id}->[0] , $germline_stat{$id}->[3]/$germline_stat{$id}->[1]);
	}
	my ($ave_abund_r ,$ave_uniq_r ,$sub_abund_r ,$sub_uniq_r,$sub_adund_r_sum ,$sub_uniq_r_sum) = (0,0,0,0,0,0);# calculate rate
	my $uniq_num = scalar(keys %germline_stat);
	$ave_abund_r = $germline_stat{$id}->[0]/$sum_seq_a*$uniq_num*100;
	$ave_uniq_r = $germline_stat{$id}->[1]/$sum_seq_u*$uniq_num*100;
	$sub_abund_r = $germline_stat{$id}->[2]/$sum_sub_a*$uniq_num*100;
	$sub_uniq_r = $germline_stat{$id}->[3]/$sum_sub_u*$uniq_num*100;
	$sub_adund_r_sum = $germline_stat{$id}->[2]/$germline_stat{$id}->[4]*100;
	$sub_uniq_r_sum = $germline_stat{$id}->[3]/$germline_stat{$id}->[5]*100;
	my $more5_rate = 0;
	$germline_stat{$id}->[6] = 0 if(!defined($germline_stat{$id}->[6]));
	$more5_rate = $germline_stat{$id}->[6]/$germline_stat{$id}->[1]*100;
	my $ave_more5_r = $germline_stat{$id}->[6]/$sum_more5*$uniq_num*100;


	print "$id\t$germline_stat{$id}->[0]\t$germline_stat{$id}->[1]\t$germline_stat{$id}->[2]\t$germline_stat{$id}->[3]\t$abund_r\t$uniq_r\t$ave_abund_r\t$ave_uniq_r\t$sub_abund_r\t$sub_uniq_r\t$germline_stat{$id}->[6]\t$more5_rate\t$ave_more5_r\n";
	
	# filter $interval_filter ,$ave_r_filter
#	next if($germline_stat{$id}->[0] <50 || $germline_stat{$id}->[1] <50);
	if($abund_r >= $interval_filter && $uniq_r >= $interval_filter && $sub_abund_r >= $ave_r_filter && $sub_uniq_r >= $ave_r_filter && $more5_rate>=$filter_more5_1 && $ave_more5_r>=$filter_more5_2)
	{
		print O "$id\n$backup{$id}\n";
	}

}



