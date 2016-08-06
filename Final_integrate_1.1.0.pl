#!/usr/bin/perl -w
use strict;
use Data::Dumper;

die "perl $0 <in.fa><out.fa><gene><mis_num><annontion.file><final.out.fa><reference.file>\n" unless(@ARGV==7);

# print out format:
# # >V_8:1:11:1435:1357:1
# # V_seed_name:flag:combined_seq_num:abundance:unique:major(1) or minor(2)
#

my $gene = $ARGV[2];
my $mis_num = $ARGV[3];
my $anno_file = $ARGV[4];
my $final_out = $ARGV[5];
my $ref_file = $ARGV[6];

open O ,">$ARGV[1]" or die;
print O "#ID\tScore\tTotal_reads_support\tUnique_reads_support\tCombined_seq\tRecommend_Gene/allele\tFlag\tMapped_to_known_germline\tMapped_identity\tMapped_mismatch\tMapped_Deviated_bases\n";

open F, ">$final_out" or die;
#my $major_rate = 10;
my $mini_edge = 5;
my $mask_dege = 0;

#--------	read reference file	-----------------
my %known_gene_name;
open R, "$ref_file" or die;
$/=">";
<R>;
while(<R>)
{
	chomp;
	my @line = split /\s+/,$_;
	my ($id , $species) = (split /_/,$line[0])[0,-1];
	my ($id_gene, $id_allele) = split /\*/,$id;
	$id_allele =~ s/^0//;
	$known_gene_name{$species}{$id_gene}{$id_allele} = 1;
}
close R;
$/="\n";

#print Dumper(\%known_gene_name);

#-0------------		read annotation file from blast result	-------------
my %annotation;
open A, "$anno_file" or die;
<A>;
while(<A>)
{
	chomp;
	my @line = split;

	my $tag;
	my $naming;

	if($line[7]==0){
		$tag = "known";
		$naming = $line[-1];
	}elsif(($gene eq "V" && $line[7]<=7) || ($gene eq "J")){
		$tag = "novel_allele";
		my ($id , $species) = (split /_/,$line[-1])[0,-1];
		my ($id_gene, $id_allele) = split /\*/,$id;
		my $max = (sort {$b <=> $a} keys %{$known_gene_name{$species}{$id_gene}})[0];
		$max++;
#		$naming = $id_gene."*".($max);
		$known_gene_name{$species}{$id_gene}{$max} = 1;
		$max = "0$max" if(length($max)==1);
		$naming = $id_gene."*".($max);
	}else{
		$tag = "novel_gene";
	}
	
	$annotation{$line[3]} = [("$tag\t$line[-1]\t$line[0]\t$line[7]\t$line[1]",$tag,$naming)];
}
close A;


#-----	read fa file	----------
my %all;
open I, "$ARGV[0]" or die;
while(<I>)
{
	chomp;
	s/>//;
	my $id = $_;
	my ($fir,$abund,$tot_s,$uniq_s,$comb) = (split /:/,$id)[0,1,2,3,4];
	chomp(my $seq = <I>);
	if($annotation{$id}->[1] eq "known" || $annotation{$id}->[1] eq "novel_allele"){
		print O "$fir\t$abund\t$tot_s\t$uniq_s\t$comb\t$annotation{$id}->[2]\t$annotation{$id}->[0]\n";
	}else{
		$all{$id} = [($seq,$abund)];
	}

#	$all{$id} = [($seq,$abund)];
	print F ">$fir\n$seq\n";
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

my $naming_start = 65;
#------	    filter	---------
#open O ,">$ARGV[1]" or die;
#print O "#ID\tScore\tTotal_reads_support\tUnique_reads_support\tCombined_seq\tRecommend_Gene/allele\tFlag\tMapped_to_known_germline\tMapped_identity\tMapped_mismatch\tMapped_Deviated_bases\n";

while(1)
{
	last if(scalar keys %matrix == 0);

	my @left_ids = keys %matrix;
	my $id1 = (sort {$all{$b}->[1]<=>$all{$a}->[1]} @left_ids)[0];

#	my $id1 = (sort {scalar @{$matrix{$b}}<=>scalar @{$matrix{$a}}} keys %matrix)[0];
#	my $max_abund = 0;
#	for my $id2(@{$matrix{$id1}}){
#		$max_abund = $all{$id2}->[1] if($max_abund < $all{$id2}->[1]);
#	}
	
#	my %new_t;
	my %del_id;
#	$del_id{$id1} = 1;
#
#	#-----------	print out	-----------------
#
	my @new = sort {$all{$b}->[1]<=>$all{$a}->[1]} @{$matrix{$id1}};
#	unshift @new, $id1;
	unshift @new, "NA";

	for(my $i= 1 ; $i<=$#new ; $i++){
		my @detail = split /:/,$new[$i];
		my $recom_gene = "$gene-".chr($naming_start)."*0$i";

		print O "$detail[0]\t$detail[1]\t$detail[2]\t$detail[3]\t$detail[4]\t$recom_gene\t$annotation{$new[$i]}->[0]\n";

#		print F ">$detail[0]\n$all{$new[$i]}->[0]\n";
		$del_id{$new[$i]} = 1;
	}
	$naming_start++;
#	for my $id2(@{$matrix{$id1}}){
#		$del_id{$id2} = 1;
#		my $rate = $max_abund/$all{$id2}->[1];
#		next if($rate > $major_rate);
#		$new_t{$id2} = $rate;
#	}
#	my @new = sort{$new_t{$a}<=>$new_t{$b}} keys %new_t;

#	if(scalar @new >1){
#		print O "$new[0]:1\n$all{$new[0]}->[0]\n";
#		print O "$new[1]:2\n$all{$new[1]}->[0]\n";
#		$del_id{$new[0]} = 1;
#		$del_id{$new[1]} = 1;
#	}else{
#		print O "$new[0]:1\n$all{$new[0]}->[0]\n";
#		$del_id{$new[0]} = 1;
#	}
	
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

close O;
close F;

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


