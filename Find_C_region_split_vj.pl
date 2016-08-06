#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

unless (@ARGV == 9)
{
	print "Usage:\n";
	print "perl $0 <*.primer><*.fa><v.fa.gz><j.fa.gz><mismatch_num><miss_v><j_retain_num><num_for_C><V_seed_len>\n";
	print "primer file: sequence with forward starnd\n";
	exit;
}
my ($in_file1 , $in_file2 , $out1 , $out2 , $mismatch_num ,$miss_v ,$retain_j_num , $num_for_c, $v_seed_len) = @ARGV;
open IN1 , "$in_file1" or die;
if($in_file2 =~ /\.gz/){
	open IN2 , "gzip -dc $in_file2|" or die;
}
else{
	open IN2 , "$in_file2" or die;
}
open V , "|gzip>$out1" or die;
open J , "|gzip>$out2" or die;



my $indent_num = 100;

# generation base
my %degeneration = ("A"=>"A","C"=>"C","T"=>"T","G"=>"G",
		"M"=>"AC","R"=>"AG","W"=>"AT","S"=>"CG","Y"=>"CT","K"=>"GT",
		"V"=>"ACG","H"=>"ACT","D"=>"AGT","B"=>"CGT","N"=>"ACGT");


#---------------------------------------------------
# construct primers with:
# 1. mismatch<=2
# 2. delete 1 base of 5' primer


my %Primer;
while(<IN1>)
{
	chomp;
	s/^>//;
	my $id = $_;
	chomp(my $seq = <IN1>);
	if(length($seq)>=$num_for_c){
		$seq=substr($seq,0,$num_for_c);
	}

	$seq =~ tr/acgt/ACGT/;
	$seq =~ tr/ACGT/TGCA/;
	$seq = reverse $seq;

	if($seq =~ /[^ACGT]/)# has degenerate base
	{
		my %mulit;
		my @temp = split // , $seq;
		for my $base(@temp)
		{
			if(scalar keys %mulit == 0){
				$mulit{$_} = 1 for (split //,$degeneration{$base});
			}
			else{
				for my $part_seq (keys %mulit){
					$mulit{"$part_seq$_"} = 1 for (split //,$degeneration{$base});
					delete $mulit{$part_seq};
				}
			}
		}
		for (keys %mulit){
			&primer_mismatch_2($_);
			my $d_1_s = substr($_,1);
			$Primer{length($d_1_s)}{$d_1_s} = 1;
		}
	}
	else
	{
		&primer_mismatch_2($seq);
		my $d_1_s = substr($seq , 1);
		$Primer{length($d_1_s)}{$d_1_s} = 1;
	}
	
}
close IN1;

#print Dumper(\%Primer);
# create primer sequence with mismatch <=2bp
sub primer_mismatch_2
{
        my $s = shift @_;
        $Primer{length($s)}{$s} = 1;
        my %new;
        $new{$s} = 1;

        for(my $i=1 ; $i<=$mismatch_num ; $i++)
        {
                for my $s_new(keys %new)
                {
                        for(my $j=0 ; $j<length($s_new) ; $j++)
                        {
                                for("A","C","T","G"){
                                        my $s_new_1 = $s_new;
                                        substr($s_new_1,$j,1) = $_;
                                        $Primer{length($s_new_1)}{$s_new_1} = 1;
                                        $new{$s_new_1} = 1;
                                }
                        }
                }
        }

}





#-----------------------------------------
my %V_all;
my %J_all;

my ($find_p_n , $all_n) = (0,0);
my ($len_f_v,$len_f_j) = (0,0);
while(<IN2>)
{
	chomp(my $id = $_);
	chomp(my $seq = <IN2>);
#	chomp(my $l3 = <IN2>);
#	chomp(my $l4 = <IN2>);
	$all_n++;
	
	my $flag = 0;

	$indent_num = length($seq)-1 if($indent_num>=length($seq));

	for(my $i=0 ; $i<$indent_num ; $i++)
	{
		for my $len(keys %Primer){
			my $l =substr($seq , $i , $len);
			next if(length($l)<$len);
			if(exists $Primer{$len}{substr($seq , $i , $len)}){
				$find_p_n++;
				$flag = 1;
				$seq = substr($seq , $i+$len);
				$seq =~ tr/ACGT/TGCA/;
				$seq = reverse $seq;
				my $s_v = substr($seq , 0 , length($seq)-$miss_v);
				my $s_j = substr($seq , -$retain_j_num);
				if(length($s_v)>=$v_seed_len){
#					print V "$id\n$s_v\n";
					$V_all{$s_v}++;
				}
				else{$len_f_v++;}
				if(length($s_j)==$retain_j_num){
#					print J "$id\n$s_j\n";
					$J_all{$s_j}++;
				}else{$len_f_j++;}
				last;
			}
		}
		last if($flag);
	}

	next if($flag);

	my $reverse_seq = reverse $seq;
	$reverse_seq =~ tr/ACGT/TGCA/;
	$flag = 0;
	for(my $i=0 ; $i<$indent_num ; $i++)
	{
		for my $len(keys %Primer){
			my $l =substr($seq , $i , $len);
			next if(length($l)<$len);
			if(exists $Primer{$len}{substr($reverse_seq , $i , $len)}){
				$find_p_n++;
				$flag = 1;
				$seq = $reverse_seq;
				$seq = substr($seq , $i+$len);
				$seq =~ tr/ACGT/TGCA/;
				$seq = reverse $seq;
				 my $s_v = substr($seq , 0 , length($seq)-$miss_v);
				 my $s_j = substr($seq , -$retain_j_num);
				 if(length($s_v)>=$v_seed_len){
#					 print V "$id\n$s_v\n";
					 $V_all{$s_v}++;
				 }
				else{$len_f_v++;}
				if(length($s_j)==$retain_j_num){
#					print J "$id\n$s_j\n";
					$J_all{$s_j}++;
				}else{$len_f_j++;}								 
				last;
			}
		}
		last if($flag);
	}
	
	next if($flag);

}

my $flag_v = 0;
for(sort {$V_all{$b}<=>$V_all{$a}} keys %V_all)
{
	$flag_v++;
	print V ">v$flag_v:$V_all{$_}\n$_\n";
}
my $flag_j = 0;
for(sort {$J_all{$b}<=>$J_all{$a}} keys %J_all)
{
	$flag_j++;
	print J ">j$flag_j:$J_all{$_}\n$_\n";
}


print "sequence num: $all_n\n";
print "find primer seq: $find_p_n\t",$find_p_n/$all_n*100,"\n";
print "no primer seq: ",$all_n-$find_p_n,"\n";
print "Retain_seq_v: ",$find_p_n-$len_f_v,"\t",($find_p_n-$len_f_v)/$all_n*100,"\n";
print "Retain_seq_j: ",$find_p_n-$len_f_j,"\t",($find_p_n-$len_f_j)/$all_n*100,"\n";
