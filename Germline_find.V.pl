#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use Graph::Undirected;



=head1 Usage
	perl Germline_find.pl

=head1 [parameters]
	
	-i	<S> input file
	-m	<I> the read number needed for seed extension [2]
	-mr	<F> Rate of the read number needed for seed extension [0.03]
	-n	<F> the unique read number needed for seed extension [2]
	-nr     <F> Rate of the unique read number needed for seed extension [0.05]
	-mn     <I> the read number needed for seed extension(non-CDR3 end)[2]
	-nn     <> the unique read number needed for seed extension(non-CDR3 end) [2]
	-l	<I> the minimum length for output [0]
	-p	<S> the prefix of id name

	print out format:
	>V_100:1:7:7
	V_seed_name:flag:abundance:unique
=cut

my ($infile , $min_num , $uniq_num, $germline_len , $name,$min_num_r , $uniq_num_r, $min_num_non_cdr3, $uniq_num_non_cdr3);

GetOptions(
		"i:s"=>\$infile,
		"m:i"=>\$min_num,
		"mr:f"=>\$min_num_r,
		"n:i"=>\$uniq_num,
		"nr:f"=>\$uniq_num_r,
		"mn:f"=>\$min_num_non_cdr3,
		"nn:f"=>\$uniq_num_non_cdr3,
		"l:i"=>\$germline_len,
		"p:s"=>\$name
	  );


die `pod2text $0` if(!$infile && !$name);

$min_num=2 unless(defined($min_num));
$min_num_r=0.03 unless(defined($min_num_r));
$uniq_num=2 unless(defined($uniq_num));
$uniq_num_r=0.05 unless(defined($uniq_num_r));
$germline_len = 0 unless(defined($germline_len));
$min_num_non_cdr3=2 unless(defined($min_num_non_cdr3));
$uniq_num_non_cdr3=2 unless(defined($uniq_num_non_cdr3));


my $raw_noncdr3_mini_mum = $min_num_non_cdr3;
#------		1. read input file	---------
# 1. store base in a hash, before the seed and after the seed
# 2. count the base number for each postion

if($infile =~/\.gz/){
	open I, "gzip -dc $infile|" or die;
}else{
	open I, "$infile" or die;
}

my $seed = <I>;
my ($abund_all ,$uniq_all);

($seed, $abund_all ,$uniq_all) = split /:/,$seed;# seed sequence
if($abund_all < $min_num or $uniq_all<$uniq_num){
	exit;
}

my %back; # store raw sequence for further analysis
my %freq; # store the base for each position
my $uniq_seq = 0;
while(<I>)
{
	chomp;
	my ($abund, $start , $end) = (split /[:-]/,$_)[1,2,3];

	chomp(my $seq = <I>);
	$back{$seq} = [($abund,$start , $end)];
	$uniq_seq++;

	my $before = substr($seq , 0 , $start-1);
	$before = reverse $before;
	my @b = split //,$before;
	my $i = 0;
	for(@b){
		$i--;
		$freq{$i}{$_} += $abund;
	}

	my $after = substr($seq , $end);
	my @a = split // , $after;
	$i=0;
	for(@a){
		$i++;
		$freq{$i}{$_} += $abund;
	}
	$freq{0}{$seed} += $abund;
}
close I;
#print Dumper(\%freq);

#-------	2. filter the bases	------------
# 1. find the start and end position(satisfy the base number required)
# 2. remove the base without satisfied

my ($min , $max) = (sort {$a<=>$b} keys %freq)[0,-1];

for my $pos(sort {$a<=>$b} keys %freq)
{
	my $flag = 0;
	for (keys %{$freq{$pos}}){
		if($pos>0 && $freq{$pos}{$_} >= $min_num_r*$abund_all && $freq{$pos}{$_} >= $min_num){
			$flag = 1;
		}elsif($pos<0 && $freq{$pos}{$_} >= $raw_noncdr3_mini_mum){
			$flag = 1;
		}else{
			delete $freq{$pos}{$_};
		}
	}
	if($flag==0 && $pos <0){
		$min = $pos+1;
		delete $freq{$pos};
	}
	elsif($flag == 0 && $pos >0){
		$max = $pos-1;
		delete $freq{$pos};
	}
}
#print "$min\t$max";
#print Dumper(\%freq);


#-------------		3. seed extension	------------
# base by base extension

my %connect;
my %connect_1;


# extension for the right of seed
$connect{$seed} = [($uniq_seq,$freq{0}{$seed})];

for(my $i=1 ; $i<=$max ; $i++)
{
        my $flag_1 = 0;
        for my $raw (keys %connect)
        {
                my $flag = 0;
                my @raw_new;
                ($flag, @raw_new) = &germline_judge($raw,$i,keys %{$freq{$i}});
                if($flag){
                        for(@raw_new){
                                $connect{$_} = [(0,0)];
                        }
                        $flag_1 = 1;
                }

                $connect_1{$raw} = [(0,0)] if($flag==0);
                delete $connect{$raw};
        }
        last if($flag_1==0);

}
%connect_1 = (%connect_1,%connect);
%connect=();

# extension for the left of seed
my %connect_2;
for(my $i=-1 ; $i>=$min ; $i--)
{
	my $flag_1 = 0;
	for my $raw (keys %connect_1)
	{
		my $flag = 0;
		my @raw_new;
		
		($flag, @raw_new) = &germline_judge($raw,$i,keys %{$freq{$i}});
		if($flag){
			for(@raw_new){
				$connect_1{$_} = [(0,0)];
			}
			$flag_1 = 1;
		}

		$connect_2{$raw} = [(0,0)] if($flag==0);
		delete $connect_1{$raw};

	}
	last if($flag_1==0);
}
%connect_2 = (%connect_2,%connect_1);
%connect_1=();


#print Dumper(\%connect_2);
#------------		4. output sequence	----------------

# get the final abundance and unique number
for my $new (keys %connect_2)
{
	if(length($new) >= $germline_len){
		($connect_2{$new}->[0],$connect_2{$new}->[1]) = &find_abund_uniq($new);
	}else{
		delete $connect_2{$new};
	}
}

my $f = 0;
for my $new (sort {$connect_2{$b}[1] <=> $connect_2{$a}[1]} keys %connect_2)
{
	$f++;
	print ">$name:$f:$connect_2{$new}->[1]:$connect_2{$new}->[0]\n$new\n";
}


#-----------
# judge condition:
#	surpporting  sequences' abundance and unique number
#
sub germline_judge
{
	my ($seq_for_jud,$i,@add_bases) = @_;
	
	my %add_bases_h;
	if($i<0){
		$add_bases_h{$_.$seq_for_jud} = 1 for(@add_bases);
	}else{
		$add_bases_h{$seq_for_jud.$_} = 1 for(@add_bases);
	}
		
	my $flag = 0;
	my (%uniq_n_2, %abund_2);
	my @success;
	
	my($uniq_all_new , $abund_all_new) = (0,0);
	for my $seq(keys %back)
	{
		my ($n,$s,$e) = ($back{$seq}->[0],$back{$seq}->[1],$back{$seq}->[2]);
		if($i<0){
			$s=$s+$i;
			$e=$s+length($seq_for_jud);
		}else{
			$e=$e+$i;
			$s=$e-length($seq_for_jud);
		}
		next if($s<=0 or $e >length($seq));	

		$abund_all_new+=$n;
		$uniq_all_new++;
		my $subseq = substr($seq,$s-1,$e-$s+1);
		
		
		for(keys %add_bases_h){
			if($subseq eq $_){
				$uniq_n_2{$_}++;
				$abund_2{$_} += $back{$seq}->[0];
#				if($uniq_n_2{$_}>=$uniq_num_r*$uniq_all && $uniq_n_2{$_}>=$uniq_num && $abund_2{$_} >= $min_num_r*$abund_all && $abund_2{$_} >= $min_num){
#					push @success,$_;
#					delete $add_bases_h{$_};
#					$flag = 1;
#				}
			}
		}
#		last if(scalar keys %add_bases_h == 0);

	}
        for(keys %add_bases_h){
		next if(!exists $uniq_n_2{$_});
		if($i>0 && $uniq_n_2{$_}>=$uniq_num_r*$uniq_all_new && $uniq_n_2{$_}>=$uniq_num && $abund_2{$_} >= $min_num_r*$abund_all_new && $abund_2{$_} >= $min_num){
			push @success,$_;
			$flag = 1;
		}elsif($i<0 && $uniq_n_2{$_}>=$uniq_num_r*$uniq_all_new && $uniq_n_2{$_}>=$uniq_num_non_cdr3 && $abund_2{$_} >= $min_num_r*$abund_all_new && $abund_2{$_} >= $min_num_non_cdr3){
			push @success,$_;
			$flag = 1;
		}
	}
	return($flag ,@success);
}

sub find_abund_uniq
{
	my $seq_for_jud = shift @_;
	my ($uniq_n_2, $abund_2)=(0,0);
	for my $seq(keys %back)
	{
		if($seq=~/$seq_for_jud/){
			$uniq_n_2++;
			$abund_2 += $back{$seq}->[0];
		}
	}
	return ($uniq_n_2, $abund_2);
}


