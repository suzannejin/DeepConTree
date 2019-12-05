#!/usr/bin/perl
#
################################################################################
####################
# This script read as input file the list of all datasets to be treated. It will
# therefore create
# a folder for each of them and move the corresponding fasta file in each folder
# Finally, it will create a shell script for each dataset to be used
# ARGV[0] : list of Pfam family for which a script is needed
# ARGV[1] : MSA methods to apply on the scripts
# ARGV[2] : export T-Coffee MODE
# ARGV[3] : export T-Coffee MODE $EXP (2=square; 3=cubic)
#
################################################################################
####################
use strict;
use File::Copy;

my $list=$ARGV[0];
my $msa=$ARGV[1];
my $mode=$ARGV[2];
my $exp=$ARGV[3];

open (LIST,"$list");
my $fasta;
my $template;
my $trees;
my $matrices;
my $output;

while (my $line=<LIST>)
{
        chomp $line;
        $fasta=$line.".fa";
        $template=$line."_ref.template_list2";
        $trees=$line."_unw_d4_"."$msa"."$mode"."-"."$exp.trees";
        $matrices=$line."_unw_d4_"."$msa"."$mode"."-"."$exp".".matrices";
        $output=$line."_"."$msa"."_phylo3D_unw_d4-"."$mode"."-"."$exp".".sh";
        open (SHELL,">$output");
        print SHELL "#!/bin/bash\n";
        print SHELL "cd $line\n";
        #print SHELL "export BS_SQRLEN=1\n";
        print SHELL "export THREED_TREE_MODE=$mode\n";
        print SHELL "export THREED_TREE_MODE_EXP=$exp\n";
        print SHELL "export THREED_TREE_NO_WEIGHTS=1\n";
        #print SHELL "t_coffee -other_pg seq_reformat -in $fasta -in2 $template -action +tree replicates 100 +evaluate3D distances +tree2bs first +print_replicates -output newick > $trees \n";
        print SHELL "t_coffee -other_pg seq_reformat -in $fasta -in2 $template -action +tree replicates 100 +evaluate3D distances +tree2bs first +print_replicates -output dm > $matrices \n";
        close (SHELL);
}

`chmod ug+rwx *.sh`.
close (LIST);
