# Xenacoelomorpha2019

Data from Philippe et al 2019.  TITLE. REFERENCE.

1. Sequence alignments untrimmed
  In three zip files "untrimmed_1-399", "untrimmed_400-799" and "untrimmed_800-1173".
2. Sequence alignments trimmed (individual files)
  "all_alignments_trimmed.zip"
3. Sequence alignments trimmed and concatenated from 'best' to 'worst'.
  Nexus file with character sets definining positions of each individual alignment.
  "all_genes_ordered_by_monophyly.nexus.zip"
3. Trees.
4. "monophyly.pl"  Perl script to measure support for user defined monophyletic groups.
  usage: "monophyly.pl <file with clades defined> <treefile(s) file1.tre or \*.tre>"
  requires https://metacpan.org/pod/Bio::TreeIO; https://metacpan.org/pod/Bio::Tree::Node;

Clades file has format as follows

clade_name1 TAB spp1 spp2 spp3\
clade_name2 TAB spp3 spp4 spp10

TAB between name of the clade and first spp.  Space between each spp.\
Spp can appear in more than one clade
Spp must be spelled exactly as in tree file.
