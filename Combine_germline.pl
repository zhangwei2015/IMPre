#!/usr/bin/perl -w 
use strict;
use Data::Dumper;
use File::Basename;

die "perl $0 <in> <over_min_len> <over_mat_r> <gene> <output_len_filter>\n" unless(@ARGV==5);

# print out format:
# >V_8:1:11:1435:1357
# V_seed_name:flag:combined_seq_num:abundance:unique


my $over_min_len = $ARGV[1];
my $over_mat_r = $ARGV[2];
my $edge_mask_n = 5; # the edge bases with more mismatch, V for 5' and J for 3'
my $edge_mask_m = 2;
my $gene = $ARGV[3];
my $length_restrict = $ARGV[4];

my $indir = dirname($ARGV[0]);
open B, ">$indir/$ARGV[3]_combined.log" or die;

#---	1. read the raw file	--------------------
my %all;

open I, "$ARGV[0]" or die;
while(<I>)
{
	chomp;
	s/^>//;
	my @id = split /:/,$_;
	chomp(my $seq = <I>);
	$all{"$id[0]:$id[1]"}= [($id[2],$id[3],$seq,length($seq),1)];
}
close I;


#2----	2. start to combine	-------------
my %germline = %all;
my %backup;
#print Dumper(\%germline);
my $flag_1; 
my $flag_2 = 0;
while(1)
{
	%all = ();
	%all = %germline;
	$flag_1 = scalar keys %all;# the sequence number of before combined
	%germline = ();


	for my $id1 (sort {$all{$b}->[0]<=>$all{$a}->[0]} keys %all)
	{
		next if(!exists $all{$id1});

#		print B "$id1";
		$germline{$id1} = $all{$id1}; 
		my @t = split //,$all{$id1}->[2];
		my %g1;
		for(my $i=1;$i<=scalar @t;$i++){
			$g1{$i} = [($t[$i-1],$all{$id1}->[0],$all{$id1}->[1])];
		}
		my $sam_num = 0;
		for my $id2 (keys %all)
		{
			next if($id1 eq $id2);
			my @tt = split //,$all{$id2}->[2];
			my %g2;
			for(my $i=1;$i<=scalar @tt;$i++){
				$g2{$i} = [($tt[$i-1],$all{$id2}->[0],$all{$id2}->[1])];
			}

			my ($flag,$seq,$abund,$uniq) = &merge(\%g1, \%g2);# parameter use reference, the hash will be changed
			if($flag){
#				$germline{$id1}->[0] = $abund;
#				$germline{$id1}->[1] = $uniq;
#				$germline{$id1}->[2] = $seq;
#				$germline{$id1}->[3] = length($seq);
				$sam_num+=$all{$id2}->[4];
				delete $all{$id2};
				print B "$id1\t$id2\n";
				$backup{$id1} = 1;
			}
			
		}
		
#		print B "$id1\n" if($sam_num == 0);
		delete $all{$id1};
		my ($s,$ab,$ui) = ("",0,0);
		for(sort {$a<=>$b} keys %g1){
			$s .= $g1{$_}->[0];
			$ab += $g1{$_}->[1];
			$ui += $g1{$_}->[2];
		}
		$germline{$id1}->[0] = int $ab/length($s);
		$germline{$id1}->[1] = int $ui/length($s);
		$germline{$id1}->[2] = $s;
		$germline{$id1}->[3] = length($s);
		$germline{$id1}->[4] += $sam_num;

	}
	$flag_2 = scalar keys %germline;# the sequence number of after combined

	# the combine has finished
        if($flag_1 == $flag_2)
	{
		last;
	}

}

%all=();




for (sort {$germline{$b}->[0]<=>$germline{$a}->[0]} keys %germline)
{
	if(length($germline{$_}->[2]) >= $length_restrict)# output sequence length filter
	{
		print ">$_:$germline{$_}->[4]:$germline{$_}->[0]:$germline{$_}->[1]\n$germline{$_}->[2]\n";
	}

	if(!exists $backup{$_}){
		print B "$_\n";
	}
}

close B;
sub merge
{
	my ($s1 , $s2) = @_;
	my ($seq , $uniq , $abund) = ("",0,0);
	my $flag = 0;

	my $i=1;
	my $j=scalar(keys %$s2)-$over_min_len+1;
	for(; $i<= scalar(keys %$s1)-$over_min_len+1 ; $i++)
	{
	
			my ($mat,$len,$edge_mis_n) = (0,0,0);
			# get the alignment length
			$i=1 if($j>=1);
			if($j>1){
				my $l1 = scalar(keys %$s2)-$j+1;
				my $l2 = scalar(keys %$s1);
				$len = $l1>$l2?$l2:$l1;
			}elsif($j<=1 && $i<=scalar(keys %$s1)-scalar(keys %$s2)+1){
				$j = 1;
				$len = scalar(keys %$s2)
			}else{
				$j = 1;
				$len = scalar(keys %$s1)-$i+1;
			}
			# start to align
			
			if($gene eq "V"){
				my $k=0;
				for(; $k<$len-$edge_mask_n ; $k++)
				{
					$mat++ if($$s1{$i+$k}->[0] eq $$s2{$j+$k}->[0]);
				}
				for(; $k<$len ; $k++)
				{
					$edge_mis_n++ if($$s1{$i+$k}->[0] ne $$s2{$j+$k}->[0]);
				}
			}elsif($gene eq "J"){
				my $k=0;
				for(; $k<$edge_mask_n ; $k++){
					$edge_mis_n++ if($$s1{$i+$k}->[0] ne $$s2{$j+$k}->[0]);
				}
				for(; $k<$len ; $k++)
				{
					$mat++ if($$s1{$i+$k}->[0] eq $$s2{$j+$k}->[0]);
				}
			}
			# judge whether the demand is statisfied
			
			if($mat/($len-$edge_mask_n)>=$over_mat_r && $edge_mis_n<=$edge_mask_m)
			{
				# align again, to get the base and abundance for each position
				for(my $k=0 ; $k<$len ; $k++)
				{
					if($$s1{$i+$k}->[0] eq $$s2{$j+$k}->[0]){
						$$s1{$i+$k}->[1] += $$s2{$j+$k}->[1];
						$$s1{$i+$k}->[2] += $$s2{$j+$k}->[2];
					}else{
						if($gene eq "V" && $j+$len-1< scalar keys %$s2)# $s2 is longer than $s1 at 3' end
						{
							$$s1{$i+$k} = $$s2{$j+$k};
						}elsif($gene eq "V" && $j+$len-1 == scalar keys %$s2 && $i+$len-1 == scalar keys %$s1){
							if($$s2{$j+$k}->[2]>$$s1{$i+$k}->[2]){
								$$s1{$i+$k} = $$s2{$j+$k};
							}
						}elsif($gene eq "J" && $j>1)# $s2 is longer than $s1 at 5' end
						{
							$$s1{$i+$k} = $$s2{$j+$k};
						}elsif($gene eq "J" && $j==1 && $i==1){
							if($$s2{$j+$k}->[2]>$$s1{$i+$k}->[2]){
								$$s1{$i+$k} = $$s2{$j+$k};
							}
						}

					}
				}
				
				# add the two edge bases
				if($j>1)# the $s2 is longer than $s1
				{
					my $l=0;
					for(my $k=$j-1 ; $k>=1 ; $k--){
						$$s1{$l} = $$s2{$k};
						$l--;
					}
					my($max , $min) = (sort {$b<=>$a} keys %$s1)[0,-1];
					for(my $k=$max ; $k>=$min ; $k--)
					{
						$$s1{$k+$j-1} = $$s1{$k};
						delete $$s1{$k} if($k<=0);
					}
				}elsif($j<=1 && $i>scalar(keys %$s1)-scalar(keys %$s2)+1)
				{
					
					my $l=scalar(keys %$s1)+1;
					for(my $k=$j+$len ; $k<=scalar(keys %$s2) ; $k++){
						$$s1{$l} = $$s2{$k};
						$l++;
					}
				}
				
				$flag = 1;
				# get the combined sequence, average abundance, average unique sequence
				for(sort {$a<=>$b} keys %$s1){
					$seq .= $$s1{$_}->[0];
					$abund += $$s1{$_}->[1];
					$uniq += $$s1{$_}->[2];
				}
				
				return($flag , $seq , int($abund/scalar(keys %$s1)),int($uniq/scalar(keys %$s1)));
			}
			$j--;
	}
	return(0)

}

