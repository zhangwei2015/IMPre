#!/usr/bin/perl -w
use strict;
use Data::Dumper;

die "perl $0 <in.fa><out.fa><gene>\n" unless(@ARGV==4);

# print out format:
# # >V_8:1:11:1435:1357:1
# # V_seed_name:flag:combined_seq_num:abundance:unique:major(1) or minor(2)
#

my $major_rate = 10;
my $mis_num = $ARGV[3];
my $mini_edge = 5;
my $mask_dege = 3;


my $gene = $ARGV[2];
#-----	read fa file	----------
my %all;
open I, "$ARGV[0]" or die;
while(<I>)
{
	chomp;
	my $id = $_;
	my $abund = (split /:/,$id)[-2];
	chomp(my $seq = <I>);
	$all{$id} = [($seq,$abund)];
}
close I;


# ------	cal the similar matrix	------
my %matrix;

for my $id1(keys %all)
{
	for my $id2(keys %all)
	{
#		next if($id1 eq $id2);
		my ($s1,$s2);
		if($gene eq "J"){
			$s1 = substr($all{$id1}->[0],$mask_dege);
			$s2 = substr($all{$id2}->[0],$mask_dege);
		}else{
			$s1 = substr($all{$id1}->[0],0,length($all{$id1}->[0])-$mask_dege);
			$s2 = substr($all{$id2}->[0],0,length($all{$id2}->[0])-$mask_dege);
		}
		my $flag = &check($s1,$s2);# compare two sequence
		if($flag){
			push @{$matrix{$id1}},$id2;
		}
	}
}

#print Dumper(\%matrix);


#------	    filter	---------
open O ,">$ARGV[1]" or die;
while(1)
{
	last if(scalar keys %matrix == 0);
	my $id1 = (sort {scalar @{$matrix{$b}}<=>scalar @{$matrix{$a}}} keys %matrix)[0];
	my $max_abund = 0;
	for my $id2(@{$matrix{$id1}}){
		$max_abund = $all{$id2}->[1] if($max_abund < $all{$id2}->[1]);
	}
	
	my %new_t;
	my %del_id;
	$del_id{$id1} = 1;
	for my $id2(@{$matrix{$id1}}){
		$del_id{$id2} = 1;
		my $rate = $max_abund/$all{$id2}->[1];
		next if($rate > $major_rate);
		$new_t{$id2} = $rate;
	}
	my @new = sort{$new_t{$a}<=>$new_t{$b}} keys %new_t;

	if(scalar @new >1){
		print O "$new[0]:1\n$all{$new[0]}->[0]\n";
		print O "$new[1]:2\n$all{$new[1]}->[0]\n";
		$del_id{$new[0]} = 1;
		$del_id{$new[1]} = 1;
	}else{
		print O "$new[0]:1\n$all{$new[0]}->[0]\n";
		$del_id{$new[0]} = 1;
	}
	
	# delete the sequences of @{$matrix{$id1}} in the whole %matrix
	for my $id2(@{$matrix{$id1}}){
		delete $matrix{$id2};
	}
	for my $id2(keys %matrix)
	{
		if(exists $del_id{$id2}){
			delete $matrix{$id2};
		}else{
			my @renew;
			for(@{$matrix{$id2}}){
				push @renew , $_ unless(exists $del_id{$_});
			}
			if(scalar @renew >=1){
				$matrix{$id2} = \@renew;
			}else{
				delete $matrix{$id2};
			}
		}
	}
	delete $matrix{$id1};
}



#########-------#
#		#
# check two seq #
#######

sub check
{
	my ($s1 , $s2) = @_;
	($s1 , $s2) = ($s2 , $s1) if(length($s1)<length($s2));
	
	my $j = $mini_edge ;
	for(my $i=0 ; $i<=length($s1)+$mini_edge-length($s2); $i++)
	{	
		$i = 0 if($j>=0);
		$j = 0 if($j<0);
		my $len;
		if($j>0){
			$len = length($s2)-$mini_edge;
		}elsif($j==0 && length($s2)<=length($s1)-$i){
			$len = length($s2);
		}elsif($j==0 && length($s2)>length($s1)-$i){
			$len = length($s1)-$i; 
		}

		my $sub_s1 = substr($s1,$i,$len);
		my $sub_s2 = substr($s2,$j,$len);
		my @s_1 = split //,$sub_s1;
		my @s_2 = split //,$sub_s2;

		my $err_num = 0;
		for(my $k=0 ; $k<=$#s_1 ; $k++){
			$err_num++ if($s_1[$k] ne $s_2[$k]);
		}
		if($err_num <= $mis_num){
			return 1;
		}
		$j--;
	}
	return 0;
}


