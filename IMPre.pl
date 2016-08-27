#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);


=head1	Usage

=head1 [parameters]

	-i      <S> input fa file, *.fa or *.fa.gz
	-o	<S> out directory
	-n 	<S> sample name
	
		****	data_processing:find C region	   *********
	-p      <S> input c region file, with fa format
        -c      <I> the number of necloetide acids of C region ot consider [18]
	-cm     <I> the mismatch number for C region [2]

		****	data_processing:split V and J part	**************	
	-vsplit	<I> split V sequence into multiple parts for clustering [2]
	-vm	<I> miss base number at 3' sequence for V gene [40](recommendation:TRB=40,IGH=50)
	-jm	<I> retain base number at 3 sequence for J gene [60](recommendation:TRB=60,IGH=65)

		****	clustering	**********
	-v_seed	<I> the length of seed for V clustering [200]
	-j_seed <I> the length of seed for J clustering [40]
	-v_seed_m <I> the last bases masked for V seed [20](recommendation:TRB/IGH=20,TRA/IGk/L=10)
	-j_seed_m <I> the first bases masked for J seed [5]
	-vn	<I> the number of cluster for V gene [200](recommendation:TRB=200,IGH=300)
	-jn	<I> the number of cluster for J gene [30]
	
		****	assembly:seed extension	   *************
	-v_mr	    Rate of the read number needed for V seed extension[0.15]
	-v_nr	    Rate of the unique read number needed for V seed extension[0.12]
	-j_mr       Rate of the read number needed for J seed extension [0.12]
	-j_nr	    Rate of the unique read number needed for J seed extension [0.1]
	
		****    optimization:combine && filter     *************
	-vf_ratio    <I> for V gene, ratio(#(5bp trimmed)supporing seq/#supproing seq) for filtering [1.5]
	-vf_ave2     <I> for V gene, rate2(#more 5bp/average*100) for filtering [5]
	-jf_ratio    <I> for J gene, ratio(#(5bp trimmed)supporing seq/#supproing seq) for filtering [1.5]
	-jf_ave2     <I> for J gene, rate2(#more 5bp/average*100) for filtering [5]
	
	-vf_ave      <I> for V gene, depth rate1(#(5bp trimmed)supporting/average*100) for filtering [2]
	-jf_ave      <I> for J gene, depth rate1(#(5bp trimmed)supporting/average*100) for filtering [0.5](recommendation:TRB=0.5,IGH=2)
	-v_min_e     <I> for V gene, the mismatch number for last clustering to determine marjor/minor alleles [3]
	-j_min_e     <I> for J gene, the mismatch number for last clustering to determine marjor/minor alleles [3]

	-v_lf        <I> the minimum length for output [50]
	-j_lf        <I> the minimum length for output [40]

	       ****    Annotation     *************
	-known	<S> Known germline sequences file. uesd for annotation. FASTA format,(the id should be like this: germline_name_flag_specie, such as TRBV1-1*01_F_Human)
	-re_v	<I> To recommend a gene name, V gene mismatch number allowed for clustering [7]
	-re_j	<I> To recommend a gene name, J gene mismatch number allowed for clustering [5]

=head1 Note
	
	1. If the sequence including C region, the compulsory parameters: -i -o -n -p
	2. If the sequence without C region, the sequence must be forward strand and the compulsory parameters: -i -o -n
	3. -known: inferred germline would be aligned to the provided known germline sequences. Actually, IMPre will align all inferred germline to Human and Mouse known germline sequences with default.So "-known+Human+Mouse" will be used for annotation
	4. There are lots of parameters, however, almost of them are not need to reset, just using the values by recommendation.

	Version: IMPre-1.1.0
	update: 2016.8.1
	
=cut

my ($in , $out , $pf,$vsplit_num,$vseed_num,$jseed_num,$name , $miss_v_num, $retain_j_num,$num_for_c,$mismatch_for_c, $v_seed_len,$j_seed_len,$j_reads_sup_r, $j_reads_sup_n, $j_min_out,$v_reads_sup_r, $v_reads_sup_n , $v_min_out, $vcom_ove_len, $vcom_ove_mat, $jcom_ove_len, $jcom_ove_mat);

my ($vfil_len_stat , $vf_ratio , $vf_ave, $jfil_len_stat , $jf_ratio , $jf_ave, $vseed_mask, $jseed_mask,$igh_flag);
my ($known_germ_f, $re_v_mis, $re_j_mis);

$igh_flag = 0;

my ($vf_ratio2 , $vf_ave2 , $jf_ratio2 , $jf_ave2);
my ($v_min_e,$j_min_e);

GetOptions(
		"i=s" => \$in,
		"p:s" => \$pf,
		"o=s" => \$out,
		"n=s" => \$name,
		"vsplit:i" => \$vsplit_num,
		"vn:i" => \$vseed_num,
		"jn:i" => \$jseed_num,
		"vm:i" => \$miss_v_num,
		"jm:i" => \$retain_j_num,
		"c:i" => \$num_for_c,
		"cm:i" => \$mismatch_for_c,
		"v_seed:i" => \$v_seed_len,
		"j_seed:i" => \$j_seed_len,
		"j_mr:f" => \$j_reads_sup_r,
		"j_nr:f" => \$j_reads_sup_n,
		"j_lf:i" => \$j_min_out,
		"v_mr:f" => \$v_reads_sup_r,
		"v_nr:f" => \$v_reads_sup_n,
		"v_lf:i" => \$v_min_out,
		"vfil_l:i" => \$vfil_len_stat,
		"vf_ratio:f" => \$vf_ratio,
		"vf_ave:f" => \$vf_ave,
		"jfil_l:i" => \$jfil_len_stat,
		"jf_ratio:f" => \$jf_ratio,
		"jf_ave:f" => \$jf_ave,
		"-v_seed_m:i" =>\$vseed_mask,
		"-j_seed_m:i" =>\$jseed_mask,
		"vf_ave2:f" => \$vf_ave2,
		"jf_ave2:f" => \$jf_ave2,
		"v_min_e:i" => \$v_min_e,
		"j_min_e:i" => \$j_min_e,
		"known:s" => \$known_germ_f,
		"re_v:i" => \$re_v_mis,
		"re_j:i" => \$re_j_mis
	  );


die `pod2text $0` if (!$in or !$out  or !$name);

$vsplit_num = 2 unless(defined($vsplit_num));
$vseed_num = 200 unless(defined($vseed_num));
$jseed_num = 30 unless(defined($jseed_num));
$miss_v_num = 40 unless(defined($miss_v_num));
$retain_j_num = 60 unless(defined($retain_j_num));
$num_for_c = 18 unless(defined($num_for_c));
$mismatch_for_c = 2 unless(defined($mismatch_for_c));


$v_seed_len = 200 unless(defined($v_seed_len));
$vseed_mask = 20 unless(defined($vseed_mask));
$j_seed_len = 40 unless(defined($j_seed_len));
$jseed_mask = 5 unless(defined($jseed_mask));

#--- extension	-----------
$j_reads_sup_r = 0.12 unless(defined($j_reads_sup_r));
$j_reads_sup_n = 0.1 unless(defined($j_reads_sup_n));
$v_reads_sup_r = 0.15 unless(defined($v_reads_sup_r));
$v_reads_sup_n = 0.12 unless(defined($v_reads_sup_n));

#----  combine && filter  --------------
$vcom_ove_len = int($v_seed_len*0.5+10);
$jcom_ove_len = int($j_seed_len*0.5+5);
$v_min_e = 3 unless(defined($v_min_e));
$j_min_e = 3 unless(defined($j_min_e));

#$vfil_len_stat = 50 unless(defined($vfil_len_stat));
$vf_ratio = 1.5 unless(defined($vf_ratio));
$vf_ave = 2 unless(defined($vf_ave));
$vf_ave2 = 5 unless(defined($vf_ave2));

$v_min_out = 50 unless(defined($v_min_out));
#$jfil_len_stat = 40 unless(defined($jfil_len_stat));
$jf_ratio = 1.5 unless(defined($jf_ratio));
$jf_ave = 0.5 unless(defined($jf_ave));
$j_min_out = 40 unless(defined($j_min_out));
$jf_ave2 = 5 unless(defined($jf_ave2));

#------	    annotation	-----------
$re_v_mis = 7 unless(defined($re_v_mis));
$re_j_mis = 5 unless(defined($re_j_mis));


#----- change to absolute path
my $current_path=$ENV{"PWD"};
$in = "$current_path/$in" unless($in=~/^\//);
if(defined($pf)){
	$pf = "$current_path/$pf" unless($pf=~/^\//);
}
$out = "$current_path/$out" unless($out=~/^\//);


open A, ">$out/Execute_all.sh" or die;
#   1. data processing   -------------
open O, ">$out/${name}_data_processing.sh" or die;

if(defined($pf)){
	print O "perl $Bin/Find_C_region_split_vj.pl $pf $in $out/${name}_V_all.fa.gz $out/${name}_J_all.fa.gz $mismatch_for_c $miss_v_num $retain_j_num $num_for_c $v_seed_len >$out/data_processing.stat\n";
}else{
	print O "perl $Bin/Only_split_vj.pl $in $out/${name}_V_all.fa.gz $out/${name}_J_all.fa.gz $miss_v_num $retain_j_num $v_seed_len >$out/data_processing.stat\n";
}
print O "perl $Bin/split_into_part.pl $out/${name}_V_all.fa.gz $vsplit_num $out ${name}_V\n";

#if( -e "$out/V_gene_raw"){print O "rm $out/V_gene_raw/seed*.gz\n";}else{mkdir("$out/V_gene_raw",oct("0755"));}
#if( -e "$out/J_gene"){print O "rm $out/J_gene/seed*.gz\n";}else{mkdir("$out/J_gene",oct("0755"));}
#if( -e "$out/V_gene"){print O "rm $out/V_gene/seed*.gz\n";}else{mkdir("$out/V_gene",oct("0755"));}


close O;
print A "sh $out/${name}_data_processing.sh\n";
print A "echo \"$out/${name}_data_processing.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";

#	2. clustering	-------------------
open J, ">$out/${name}_J_clustering.sh" or die;
print J "f=\$(ls $out/J_gene/seed*.gz 2> /dev/null | wc -l)\nif [ \"\$f\" != \"0\" ]\nthen rm $out/J_gene/seed*.gz\nfi\nif [ ! -d \"$out/J_gene\" ]\nthen mkdir $out/J_gene\nfi\n";

my $x=$jseed_num*100;
#print J "perl $Bin/Seed_find.pl -i $out/${name}_J_all.fa.gz -l $j_seed_len -n $jseed_num  -o $out/J_gene -x $x -p 1 -g J -jm $jseed_mask\n";
print J "$Bin/Seed_ClusterV2.0 -i $out/${name}_J_all.fa.gz -l $j_seed_len -n $jseed_num  -o $out/J_gene -x $x -p 1 -g J -J $jseed_mask\n";
print J "ls $out/J_gene/seed.*.seq.gz|awk -F\".\" '{print \"mv \"\$0\" -f $out/J_gene/seed.\"\$(NF-2)\".seq.gz\";}'|sh -\n";
print J "ls $out/J_gene/seed.*.seq.gz|perl -ne 'chomp;my \@i=split /\\//,\$_;my \@t=split /\\./,\$i[-1];if(\$.==1){print \"perl $Bin/Germline_find.J.pl -i \$_ -p J_\$t[1] -mr $j_reads_sup_r -nr $j_reads_sup_n  >$out/${name}_J_Germline.all.out.raw.fa\\n\"}else{print \"perl $Bin/Germline_find.J.pl -i \$_ -p J_\$t[1] -mr $j_reads_sup_r -nr $j_reads_sup_n  >>$out/${name}_J_Germline.all.out.raw.fa\\n\"}' >$out/${name}_J_assembly.sh\n";
close J;
print A "sh $out/${name}_J_clustering.sh\n";
print A "echo \"$out/${name}_J_clustering.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";


my $n = $vseed_num*5;
$x=$vseed_num*200;
open V, ">$out/${name}_V_clustering.sh" or die;
print V "f=\$(ls $out/V_gene/seed*.gz 2> /dev/null | wc -l)\nif [ \"\$f\" != \"0\" ]\nthen rm $out/V_gene/seed*.gz\nfi\nif [ ! -d $out/V_gene ]\nthen mkdir $out/V_gene\nfi\n";
print V "f=\$(ls $out/V_gene_raw/seed*.gz 2> /dev/null | wc -l)\nif [ \"\$f\" != \"0\" ]\nthen rm $out/V_gene_raw/seed*.gz\nfi\nif [ ! -d $out/V_gene_raw ]\nthen mkdir $out/V_gene_raw\nfi\n";
for(1..$vsplit_num)
{
#	print V "perl $Bin/Seed_find.pl -i $out/${name}_V.$_.gz -l $v_seed_len -n $n -o $out/V_gene_raw -x $x -p $_  -g V -vm $vseed_mask\n";
	print V "$Bin/Seed_ClusterV2.0 -i $out/${name}_V.$_.gz -l $v_seed_len -n $n -o $out/V_gene_raw -x $x -p $_  -g V -V $vseed_mask\n";
}
print V "ls $out/V_gene_raw/*.seq.gz| perl $Bin/cat_seed.pl - $out/V_gene $vseed_num\n";
print V "ls $out/V_gene/seed.*.seq.gz|perl -ne 'chomp;my \@i=split /\\//,\$_;my \@t=split /\\./,\$i[-1];if(\$.==1){print \"perl $Bin/Germline_find.V.pl -i \$_ -p V_\$t[1] -mr $v_reads_sup_r -nr $v_reads_sup_n  >$out/${name}_V_Germline.all.out.raw.fa\\n\"}else{print \"perl $Bin/Germline_find.V.pl -i \$_ -p V_\$t[1] -mr $v_reads_sup_r -nr $v_reads_sup_n  >>$out/${name}_V_Germline.all.out.raw.fa\\n\"}' >$out/${name}_V_assembly.sh\n";

close V;
print A "sh $out/${name}_V_clustering.sh\n";
print A "echo \"$out/${name}_V_clustering.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";

#	3.1 J assembly, seed extension	---------------------

print A "sh $out/${name}_J_assembly.sh\n" ; 
print A "echo \"$out/${name}_J_assembly.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";

#	4.1 J Optimization, filtering and combination-----
open J, ">$out/${name}_J_optimization.sh" or die;

print J "perl $Bin/Filter.new.pl $out/J_gene $out/${name}_J_Germline.all.out.raw.fa 5  $out/J_Germline.all.out.fa.filter $jf_ratio 0 0 $jf_ave2 J >$out/J_filter.log\n";
print J "perl $Bin/change_id.pl $out/J_Germline.all.out.fa.filter $out/J_filter.log $out/J_Germline.all.out.fa.filter.new\n";
print J "perl $Bin/Combine_germline.pl $out/J_Germline.all.out.fa.filter.new $jcom_ove_len 1 J $j_min_out >$out/J_Germline.all.out.combine.fa\n";
print J "perl $Bin/Filter_by_depth.pl $out/J_Germline.all.out.combine.fa $jf_ave >$out/J_Germline.all.out.combine.fa.filter\n";
print J "perl $Bin/Determine_mar_min_allele.pl $out/J_Germline.all.out.combine.fa.filter $out/${name}_J_Germline.filter.fa J $j_min_e\n";
print J "perl $Bin/Score_candidate.pl $out/${name}_J_Germline.filter.fa $out/J_filter.log $out/J_combined.log $out/J_Germline.all.out.combine.fa.filter >$out/${name}_J_Germline.score.fa\n";
close J;
#print J "cat $out/J_gene/*.gz|perl $Bin/Filter.J.pl - $out/J_gene/Germline.all.out.combine.fa 5 $jfil_len_stat $out/J_gene/all_germline_out.combine.fa.filter $jf_ratio $jf_ave >$out/J_gene/filter.txt\n";
print A "sh $out/${name}_J_optimization.sh\n";
print A "echo \"$out/${name}_J_optimization.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";


#	3.2 V assembly, seed extension    ---------------------
print A "sh $out/${name}_V_assembly.sh\n";
print A "echo \"$out/${name}_V_assembly.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";


#       4.1 V Optimization, filtering and combination	------------
open V, ">$out/${name}_V_optimization.sh" or die;
print V "perl $Bin/Filter.new.pl $out/V_gene $out/${name}_V_Germline.all.out.raw.fa 5  $out/V_Germline.all.out.fa.filter $vf_ratio 0 0 $vf_ave2 V >$out/V_filter.log\n";

print V "perl $Bin/change_id.pl $out/V_Germline.all.out.fa.filter $out/V_filter.log $out/V_Germline.all.out.fa.filter.new\n";
print V "perl $Bin/Combine_germline.pl $out/V_Germline.all.out.fa.filter.new $vcom_ove_len 1 V $v_min_out >$out/V_Germline.all.out.combine.fa\n";
print V "perl $Bin/Filter_by_depth.pl $out/V_Germline.all.out.combine.fa $vf_ave >$out/V_Germline.all.out.combine.fa.filter\n";
print V "perl $Bin/Determine_mar_min_allele.pl $out/V_Germline.all.out.combine.fa.filter $out/${name}_V_Germline.filter.fa V $v_min_e\n";
print V "perl $Bin/Score_candidate.pl $out/${name}_V_Germline.filter.fa $out/V_filter.log $out/V_combined.log $out/V_Germline.all.out.combine.fa.filter >$out/${name}_V_Germline.score.fa\n";

close V;

print A "sh $out/${name}_V_optimization.sh\n";
print A "echo \"$out/${name}_V_optimization.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";

#       5. Annotation   ------------------
open O, ">$out/${name}_annotation.sh" or die;
my $know_germ_all_v = "$Bin/Annotation/Human_Mousea_TCR_BCR_F_ORF_allP_2015.07.V.fa";
my $know_germ_all_j = "$Bin/Annotation/Human_Mousea_TCR_BCR_F_ORF_allP_2015.07.J.fa";
if(defined($known_germ_f))
{
	print O "if [ ! -d $out/Known_germl ]\nthen mkdir $out/Known_germl\nfi\n";
	print O "cat $Bin/Annotation/Human_Mousea_TCR_BCR_F_ORF_allP_2015.07.V.fa $known_germ_f >$out/Known_germl/All_V_germline.fa\n$Bin/Annotation/formatdb -i $out/Known_germl/All_V_germline.fa -p F\n";
	print O "cat $Bin/Annotation/Human_Mousea_TCR_BCR_F_ORF_allP_2015.07.J.fa $known_germ_f >$out/Known_germl/All_J_germline.fa\n$Bin/Annotation/formatdb -i $out/Known_germl/All_J_germline.fa -p F\n";
	$know_germ_all_v = "$out/Known_germl/All_V_germline.fa";
	$know_germ_all_j = "$out/Known_germl/All_J_germline.fa";
}
print O "$Bin/Annotation/blastall -p blastn -d $know_germ_all_v -i $out/${name}_V_Germline.score.fa -m 8 -o $out/V_all_germline_out.fa.blast.m8\nperl $Bin/Annotation/Annotation_by_Human_Mouse.pl $out/V_all_germline_out.fa.blast.m8 $out/${name}_V_Germline.score.fa $out/${name}_V_Germline.annotation.details $know_germ_all_v\n";
print O "$Bin/Annotation/blastall -p blastn -d $know_germ_all_j -i $out/${name}_J_Germline.score.fa -m 8 -o $out/J_all_germline_out.fa.blast.m8\nperl $Bin/Annotation/Annotation_by_Human_Mouse.pl  $out/J_all_germline_out.fa.blast.m8 $out/${name}_J_Germline.score.fa $out/${name}_J_Germline.annotation.details $know_germ_all_j\n";
print O "perl $Bin/Final_integrate_1.1.0.pl $out/${name}_V_Germline.score.fa $out/${name}_V_Germline.annotation.summary V $re_v_mis $out/${name}_V_Germline.annotation.details $out/${name}_V_Germline.final.fa $know_germ_all_v\n";
print O "perl $Bin/Final_integrate_1.1.0.pl $out/${name}_J_Germline.score.fa $out/${name}_J_Germline.annotation.summary J $re_j_mis $out/${name}_J_Germline.annotation.details $out/${name}_J_Germline.final.fa $know_germ_all_j\n";
close O;

print A "sh $out/${name}_annotation.sh\n";
print A "echo \"$out/${name}_annotation.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";


#print V "cat $out/V_gene/*.gz|perl $Bin/Filter.V.pl - $out/V_gene/Germline.all.out.combine.fa 5 $vfil_len_stat $out/V_gene/all_germline_out.combine.fa.filter $vf_ratio $vf_ave >$out/V_gene/filter.txt\n";

open O, ">$out/${name}_rm_intermediate.file.sh";
print O "rm -r $out/J_gene $out/V_gene_raw $out/V_gene\n";
print O "rm $out/[VJ]_Germline.all.out.fa.filter $out/[VJ]_Germline.all.out.fa.filter.new $out/[VJ]_Germline.all.out.combine.fa\n";
print O "rm $out/${name}_[VJ]_Germline.all.out.raw.fa $out/[VJ]_Germline.all.out.combine.fa.filter $out/[VJ]_filter.log\n";
print O "rm $out/${name}_V.*.gz $out/${name}_[VJ]_all.fa.gz $out/${name}_[VJ]_Germline.filter.fa $out/${name}_[VJ]_Germline.score.fa $out/[VJ]_combined.log $out/[VJ]_all_germline_out.fa.blast.m8\n";
if(defined($known_germ_f))
{
	print O "rm -r $out/Known_germl\n";
}

close O;
print A "sh $out/${name}_rm_intermediate.file.sh\n";
print A "echo \"$out/${name}_rm_intermediate.file.sh end: \" \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\n";
close A;
