#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

unless (@ARGV == 6)
{
	print "Usage:\n";
	print "perl $0 <*.fa><v.fa.gz><j.fa.gz><miss_v><j_retain_num><V_seed_len>\n";
	print "<*.fa> must be forward strand\n";
	exit;
}
my ( $in_file2 , $out1 , $out2 , $miss_v ,$retain_j_num ,  $v_seed_len) = @ARGV;
if($in_file2 =~ /\.gz/){
	open IN2 , "gzip -dc $in_file2|" or die;
}
else{
	open IN2 , "$in_file2" or die;
}
open V , "|gzip>$out1" or die;
open J , "|gzip>$out2" or die;


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
print "Retain_seq_v: ",$all_n-$len_f_v,"\t",($all_n-$len_f_v)/$all_n*100,"\n";
print "Retain_seq_j: ",$all_n-$len_f_j,"\t",($all_n-$len_f_j)/$all_n*100,"\n";
