#!/usr/bin/perl -w

use strict;
use Bio::TreeIO;
use Bio::Tree::Node;

my @tax_array;
my $cladename;
my @array_of_clades;
my %tax_hash;
my $treefile;
my $present_in_tree;
my $count_node;
my $tree;
my $clade;
my $ingroupies;
my $in;
my $root;
my $paraphyletic;
my $distances;
my $key;
my $format;
my @nodes;
my %present_in_tree;
my $count_tree;
my $input;
my $node;
my @ingroupies;
my %ingrouphash;
my @out;
my @in;
my $sp;
my $all_spp;
my $oldroot;
my $outgroup;
my $out;
my $lca;
my $totaldistance;
my $numberofspp;
my %distancehash;
my $distance;
my $id;
my $averagedistance;
my $sqtotal;
my $std;
my %adding_up;
my $percent;
my $nodecount = 0;
my $clade_member;
my $child;
my %clade_yes;
my @sorted_nodes;
my $clademember;
my %clade_no;
my %find_bad_nodes;
my %ingroupies_in_clade;
my %to_delete;
my $aa_name;
my $aa_seq;
my @aa_seq;
my $seqfile;
my $outfile;
my $proportion;

########################################################################################################################################
#process the file that tells you which are supposed monophyletic groups
########################################################################################################################################

print "usage: monophyly.pl <file with clades defined> <treefile(s) e.g. \\*.tre>\n
clades file has format as follows

clade_name1 TAB spp1 spp2 spp3
clade_name2 TAB spp3 spp4 spp10

TAB between name of the clade and first spp.  Space between each spp.
Spp can appear in more than one clade
spp must be spelled exactly as in tree file.

treefile should be newick format
";


my $ingroups = shift;
chomp $ingroups;

open (FILE, $ingroups)  || die print "cannot open file $ingroups\n";
print "\tN_spp\t";
while (<FILE>) {
	chomp;
	if ($_ =~ m/\S/i) {
		(@tax_array) = split " ";
		$cladename = shift @tax_array;
		push @array_of_clades, $cladename;
		foreach (@tax_array) {
			push (@{$tax_hash{$cladename}}, $_);
		}
		print "$cladename\t";
		$adding_up{$cladename}{positive}=0;
		$adding_up{$cladename}{negative}=0;

		$cladename = '';
		@tax_array = ();
	}
}
print "\n";
close FILE;


########################################################################################################################################
#MAIN loop through all the treefiles
########################################################################################################################################

my $treefiles = shift;
chomp $treefiles;
my @treefiles = glob $treefiles;


foreach $treefile (@treefiles) {
#	print STDERR $treefile, "\n";
	$seqfile = $treefile;
	$seqfile =~ s/\.tre//;
	$outfile = $seqfile . '.mono';
# 	print "$seqfile\n";	
 	print "\n$treefile\t";
 	
	($present_in_tree, $count_node, $tree) = tree_process ($treefile);# $present_in_tree is ref to hash %present_in_tree
	print "$count_node";

	foreach $clade (@array_of_clades) {
		($ingroupies, $in, $root) = for_each_clade ($clade, $present_in_tree, $count_node, $tree);
	}
	$tree->cleanup_tree;	
	$input = '';
	@nodes =();
	$tree ='';
}

print "\n\n";

########################################################################################################################################
#Process tree.  Find out what spp present. 
########################################################################################################################################


sub tree_process {
	$treefile = shift;
	$format = 'newick';
	open (CHECKFORMAT, $treefile);
	while (<CHECKFORMAT>) {
		chomp;
		if ($_ =~ m/\#NEXUS/i) {
			$format = 'nexus'
		}
	}
########################################################################################################################################
#empty all the useful hashes and arrays for new tree
########################################################################################################################################

	@nodes=();
	%present_in_tree=();
	$count_tree = 0;
	$count_node = 0;

########################################################################################################################################
#set up newick/new hampshire format and nodes objects
########################################################################################################################################

	$input = new Bio::TreeIO(-file   => "$treefile",
								-format => "$format");
								
	$node = Bio::Tree::Node->new();

########################################################################################################################################
#read in next tree and its nodes
########################################################################################################################################

	while ($tree = $input->next_tree) {
		if ($count_tree > 0) {
			next;
		}
		else {
			$count_tree++;
		}
		
		@nodes = $tree->get_nodes();


		foreach $node (@nodes){
			$present_in_tree{$node->id} = 1 if ($node->is_Leaf);
			$count_node++ if ($node->is_Leaf);
		}
		return (\%present_in_tree, $count_node, $tree);
	}
}

########################################################################################################################################
#Loop through clades from clades file.
########################################################################################################################################

########################################################################################################################################
#get spp you want to be monophyletic.  Find nodes and put into @in
########################################################################################################################################

sub for_each_clade {
	$clade = shift;
	$present_in_tree = shift;
	$count_node = shift;
	$tree = shift;
	$clade_member ='';
	$child='';
	$clademember='';
	@ingroupies =();
	%ingrouphash=();
	@out=();
	@in=();
	@sorted_nodes=();
	$nodecount = 0;
	%clade_yes=();
	%clade_no=();
	%find_bad_nodes=();
	%ingroupies_in_clade = ();
	%to_delete = ();
	
	foreach (@{$tax_hash{$clade}}) {
		if ($_ =~ m/\S/i) {
			$_ =~ s/\s//;
			if ($$present_in_tree{$_}) {
				push @ingroupies, $_;
			}
		}
	}

	my $total_present = @ingroupies;

	if ($total_present <3) {
		print "\t-";
#		print "\t-\t-\t";
		return;
	}

	
	foreach $sp (@ingroupies) {
		push @in, $tree->find_node($sp);
		$ingrouphash{$sp} = 1;
	}

	foreach $all_spp (keys %present_in_tree) {
		unless (exists $ingrouphash{$all_spp}) {
			push @out, $all_spp;

		}
	}

########################################################################################################################################
#get outgroup. Find node and put into @out and reroot tree
########################################################################################################################################

	if (scalar @out > 0) {
		unless (defined $oldroot) { $oldroot = ''}; #needs to be defined for later comp
		$outgroup = $out[-1];
		chomp $outgroup;
		
		if ($outgroup ne $oldroot) { #only change if need to change
			$out = $tree->find_node($outgroup);
			$tree->reroot($out);
			$root = $tree->get_root_node;
		}
# 		else {
#  			$out = $tree->find_node($oldroot);
#  			$tree->reroot($out);
#  			$root = $tree->get_root_node;
# 		}
		$oldroot = $outgroup;
	}

	else {
		$root = '';
	}

########################################################################################################################################
#identify nodes which contain a non-ingroup member (not monophyletic)
########################################################################################################################################

	foreach $node (@nodes){
 		$clade_yes{$node->internal_id} = 0;
 		$find_bad_nodes{$node->internal_id} = 0;
		foreach $child ($node->get_all_Descendents) {
			foreach $clade_member (@out) {
				if ($child->is_Leaf and $child->id =~ m/$clade_member/) {
					$find_bad_nodes{$node->internal_id}++;
				}
			}
		}
	}

########################################################################################################################################
#for all other nodes count number of ingroup members
########################################################################################################################################
	
	foreach $node (@nodes){
		unless ($find_bad_nodes{$node->internal_id} > 0) {
			foreach $child ($node->get_all_Descendents) {
				foreach $clade_member (@ingroupies) {
					if ($child->is_Leaf and $child->id =~ m/$clade_member/) {
						$clade_yes{$node->internal_id}++;
					}
				}
			}
		}
	}

########################################################################################################################################
#sort positive nodes by how many ingroup members they have
########################################################################################################################################
	
	foreach $clademember (sort { $clade_yes{$b} <=> $clade_yes{$a} } keys %clade_yes) {
		push @sorted_nodes, $clademember;
	}


########################################################################################################################################
#Following is workaround for case when largest observed clade = 1 and there is more than one so clade with 2 or more is feasible.  
#This is not counted above as it a single member clade is not a node.  Logically, if there more than one testable members of desired ingroup 
#there must be minimum clade size of 1.
########################################################################################################################################

	if ($clade_yes{$sorted_nodes[0]} == 0 and $total_present > 1) {
		$proportion = sprintf("%.3f", 1/$total_present); 
		print "\t$proportion";# print the proportion
#		print "\t1\t$total_present\t"; #print total in biggest group THEN max size of biggest group (how many in given clade are present in tree)
	}

########################################################################################################################################
########################################################################################################################################

	else {
		$proportion = sprintf("%.3f", $clade_yes{$sorted_nodes[0]}/$total_present); 
		print "\t$proportion";# print the proportion
#		print "\t", $clade_yes{$sorted_nodes[0]}, "\t$total_present\t";#print total in biggest group THEN max size of biggest group (how many in given clade are present in tree)
	}



#	print "best node for clade $clade has ",  $clade_yes{$sorted_nodes[0]}, " out of $total_present members\n";

########################################################################################################################################
#if all the ingroup members are in largest clade it is monophyletic
########################################################################################################################################

	if ( keys %ingrouphash == $clade_yes{$sorted_nodes[0]}) {
# 		print "is monophyletic\t";
	}
	if (scalar @out == 0) {
#		print "but no outgroup members present\n";
	}
	else {
		print "";
	}

########################################################################################################################################
#get all taxa in largest monophyletic group
########################################################################################################################################

	foreach $node (@nodes) {
		if ($node->internal_id == $sorted_nodes[0]) {
			foreach $child ($node->get_all_Descendents) {
				if ($child->is_Leaf) {
					$ingroupies_in_clade{$child->id} =1 ;
				}
			}
		}
	}


	return (\@ingroupies, \@in, $root);
}

