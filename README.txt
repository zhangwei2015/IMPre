######################
Introduction

IMPre - predict BCR and TCR germline V/J genes and alleles using deep-sequencing data for rearrangement repertoire. It utilizes Clust_Seed algorithm to classify sequences and de-novel assembly for extension. V and J gene/allele are inferred seperately. FASTA format is required.

######################
System Requirement

It runs on 64-bit Linux systems.  If takes 1E6 sequences as input for example, about maximum 2.5Gb memory would be required.
Perl and C need to be installed for you system.

######################
Installation

   1. Before use it, perl(https://www.perl.org/get.html) need to be installed. 
   2. The 'Graph' module is required for the perl version.
   3. The system for C program is required.
   4. Download the IMPre to your directory, uncompress it.
   5. Make sure all programs are executable: (go to the direcotry IMPre)
   	chmod 555 *
	chmod 555 Annotation/*
	chmod 755 Test
   	

######################
Usage

   1. Create shell
   perl IMPre.pl 
	Compulsory: for sequence include C region, -i -o -n -p; for others (sequence with plus strand), -i -o -n
	Optionally: others
	all the parameters have the detail introduction if you run "perl IMPre.pl"
   this step will create multiple directory and shells

	[parameters]

            -i      <S> input fa file, *.fa or *.fa.gz
            -o      <S> out directory
            -n      <S> sample name
        
                    ****    data_processing:find C region      *********
            -p      <S> input c region file, with fa format
            -c      <I> the number of necloetide acids of C region ot consider [18]
            -cm     <I> the mismatch number for C region [2]

                    ****    data_processing:split V and J part      **************  
            -vsplit <I> split V sequence into multiple parts for clustering [4]
            -vm     <I> miss base number at 3' sequence for V gene [40](recommendation:TRB=40,IGH=50)
            -jm     <I> retain base number at 3 sequence for J gene [60](recommendation:TRB=60,IGH=65)

                    ****    clustering      **********
            -v_seed <I> the length of seed for V clustering [200]
            -j_seed <I> the length of seed for J clustering [40]
            -v_seed_m <I> the last bases masked for V seed [20](recommendation:TRB/IGH=20,TRA/IGk/L=10)
            -j_seed_m <I> the first bases masked for J seed [5]
            -vn     <I> the number of cluster for V gene [200](recommendation:TRB=200,IGH=300)
            -jn     <I> the number of cluster for J gene [30]
        
                    ****    assembly:seed extension    *************
            -v_mr       Rate of the read number needed for V seed extension[0.15]
            -v_nr       Rate of the unique read number needed for V seed extension[0.12]
            -j_mr       Rate of the read number needed for J seed extension [0.12]
            -j_nr       Rate of the unique read number needed for J seed extension [0.1]
        
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
	    -known  <S> Known germline sequences file. uesd for annotation. FASTA format(the id should be like this: germline_name_flag_specie, such as TRBV1-1*01_F_Human)
	    -re_v   <I> To recommend a gene name, V gene mismatch number allowed for clustering [7]
	    -re_j   <I> To recommend a gene name, J gene mismatch number allowed for clustering [5]

Note
            1. If the sequence including C region, the compulsory parameters: -i -o -n -p
            2. If the sequence without C region, the sequence must be forward strand and the compulsory parameters: -i -o -n
	    3. -known: inferred germline would be aligned to the provided known germline sequences. Actually, IMPre will align all inferred germline to Human and Mouse known germline sequences with default.So "-known+Human+Mouse" will be used for annotation
            4. There are lots of parameters, however, almost of them are not need to reset, just using the values by recommendation.


   2. Run shell
   2.1 it can easy to run the general sh 'Execute_all.sh': sh Execute_all.sh
   2.2 run multiple shells in seperately,so V/J could be run in parallel.
   	sh *_data_processing.sh
	sh *_J_clustering.sh
	sh *_V_clustering.sh
	sh *_J_assembly.sh
	sh *_J_optimization.sh
	sh *_V_assembly.sh
	sh *_V_optimization.sh
	sh *_annotation.sh


######################
Output
  1. output files:

	1. *_V_Germline.final.fa              : inferred germline V sequence
	2. *_V_Germline.annotation.summary    : annotation for inferred V sequence
	3. *_V_Germline.annotation.details    : details for annotation
	4. *_J_Germline.final.fa              : inferred germline J sequence
	5. *_J_Germline.annotation.summary    : annotation for inferred J sequence
	6. *_J_Germline.annotation.details    : details for annotation

  2. format for file: annotation.summary
  ID     Score   Total_reads_support     Unique_reads_support    Combined_seq Recommend_Gene/allele   Flag    Mapped_to_known_germline Mapped_identity	Mapped_mismatch Mapped_Deviated_bases
  ID: 		inferred sequence id(see *_Germline.final.fa)
  Score:	Certainty score given by IMPre, range from 0 to 100
  Total_reads_support:	the number of reads contain this inferred germline
  Unique_reads_support: the unique number of reads contain this inferred germline
  Combined_seq:	the number of extended sequences merged in the optimization step
  Recommend_Gene/allele: recommend name for the inferred germline
  Flag:		"known"(same with the known germline), "novel_allele"(<=7 mismatches for V; J:5),"novel_gene"(>7 mismathces for V; J:5)
  Mapped_to_known_germline:	nearest known germline allele
  Mapped_identity:	identity from alignment(inferred sequence && known germline sequences)
  Mapped_mismatch:	mismatch number from alignment
  Mapped_Deviated_bases: the missed bases or extra bases number at the 5' end for V( 3'end for J), comparing with nearest known germline sequence

  3. format for file: annotation.details
  Identity(%)    deviated_base   query_len       id      one_of_mapped_reference raw_identity    align_lent      mismatch        gap      query_start_position    query_end_positon       subject_start_position  subject_end_position    E-value Blast-score     Mapped_reference

######################
Testing

	directory Test/ has a data for testing
	run.sh:
	perl ../IMPre.pl -i Test.TRB.fa.gz -p TRB_C_region.txt -n T -o . -v_min_e 1 -j_min_e 1 -v_seed 40	

Please email zhangwei3@genomics.cn to report bugs/for help in installtion or usage.
