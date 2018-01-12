# CountTableProcessing
# Language: R
# Dependency: Requires plyr
# Input: prefix (for abundances and OTUs produced by Mothur)
# Output: CSV (abundance matrix)

PluMA plugin to convert OTU count output using Mothur (Schloss et al, 2009) formats 
into a CSV file of abundances.

Expected input is a file prefix, because Mothur outputs OTUs and counts in two
separate files.
The first file will be assumed to be (prefix).shared.  The .shared file in Mothur contains
raw abundances in the form of a table, where samples are rows and OTU identifiers are columns.
OTU identifiers are then mapped to their taxonomic classifications in the .taxonomy file,
which we assume to be (prefix).taxonomy:



OTU     Size    Taxonomy
Otu001  7186    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);"Porphyromonadaceae"_unclassified(100);
Otu002  4571    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);"Porphyromonadaceae"_unclassified(100);
Otu003  3897    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);"Porphyromonadaceae"_unclassified(100);



Rows of the output CSV file will then represent samples, and columns will represent OTUs.
OTU names will use their taxonomic classification, but in a simpler format.
If an OTU is classified at the genus level, only that name is used.
If an OTU is classified above the genus level, the lowest level of classification is used but with
a prefix denoting the level (K=Kingdom, P=Phylum, C=Class, O=Order, F=Family).  As an example:

Porphyromonadaceae.001 (genus Porphyromonadaceae)
F.Porphyromonadaceae.001 (family Porphyromonadaceae)

The identifier at the end distinguishes between different OTUs with the same taxonomic classification, for 
example if another OTU was classified as genus Porphyromonadaceae it would take the name Porphromonadaceae.002




