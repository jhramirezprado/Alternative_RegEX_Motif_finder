# Alternative_RegEX_Motif_finder

**Code used on this publication:**

Góngora-Castillo, E., López-Ochoa, L. A., Apolinar-Hernández, M. M., Caamal-Pech, A. M., Contreras-de la Rosa, P. A., Quiroz-Moreno, A., \***Ramírez-Prado, J. H.**, & \*O'Connor-Sánchez, A. (2020). Data mining of metagenomes to find novel enzymes: a non-computationally intensive method. *3 Biotech*, 10(2), 78.

\*corresponding author

https://doi.org/10.1007/s13205-019-2044-6

**Data mining of metagenomes to find novel enzymes: a non-computationally intensive method.**

**ABSTRACT**

Currently, there is a need of non-computationally-intensive bioinformatics tools to cope with the increase of large datasets produced by Next Generation Sequencing technologies. We present a simple and robust bioinformatics pipeline to search for novel enzymes in metagenomic sequences. The strategy is based on pattern searching using as reference conserved motifs coded as regular expressions. As a case study, we applied this scheme to search for novel proteases S8A in a publicly available metagenome. Briefly, (1) the metagenome was assembled and translated into amino acids; (2) patterns were matched using regular expressions; (3) retrieved sequences were annotated; and (4) diversity analyses were conducted. Following this pipeline, we were able to identify nine sequences containing an S8 catalytic triad, starting from a metagenome containing 9,921,136 Illumina reads. Identity of these nine sequences was confirmed by BLASTp against databases at NCBI and MEROPS. Identities ranged from 62 to 89% to their respective nearest ortholog, which belonged to phyla Proteobacteria, Actinobacteria, Planctomycetes, Bacterioidetes, and Cyanobacteria, consistent with the most abundant phyla reported for this metagenome. All these results support the idea that they all are novel S8 sequences and strongly suggest that our methodology is robust and suitable to detect novel enzymes.

**Metagenomic data**

Metagenomic sequences used in this work were obtained from the publicly available MG-RAST database [(Keegan et al. 2016](https://link.springer.com/protocol/10.1007%2F978-1-4939-3369-3_13), ID 4536384.3, and name X10-JUL09, which were obtained from mesotrophic water of a Yucatan cenote (20°90′95.67″ N/88°86′69.47″ W).

**Metagenome assembling and translation into amino acids**

Illumina raw reads were checked for quality using FastQC (v 0.11.5) (Andrews 2010) and assembled into contigs using Newbler (v2.9) or [MEGAHIT](https://github.com/voutcn/megahit) (v 1.14) [(Li et al. 2015)](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) software.

The assembled contigs were translated to amino acids (AA) into the six reading frames (+ 1, + 2, + 3, − 1, − 2, − 3) using the “transeq” tool from EMBOSS Toolkit (v 6.6.0.0) (Rice et al. 2000) by typing the following command line in a terminal window:

<pre><code>
$transeq -sequence metagenome_contigs.fna -outseq metagenome_aa_file.faa -frame 6 -Table 11 -nomethionine
</code></pre>

Before pattern searching of translated contigs, “line breaks or new lines” in the multi-fasta file were removed using the following perl command:

<pre><code>
$cat metagenome_aa_file.faa | perl -ne 'if(/^>/){print "\n",$_;next;}else{chomp;print;}' > lineal_metagenome_aa_file.faa
</code></pre>

**Pattern matching using regular expressions**

As the three motifs which shape the characteristic proteases S8A catalytic triad (CT) are well conserved, patterns for the most conserved positions were designed according to the diversity alignments reported in the conserved domains database of NCBI (https://www.ncbi.nlm.nih.gov/) and from thoroughly curated proteases S8 sequences (34 type S8A and 3 type S8B sequences) from the [MEROPS database](https://www.ebi.ac.uk/merops/) aided by the Motif discovery tool of the [MEME suite](https://meme-suite.org/meme/) (Fig. 1a) ([Bailey et al. 2015](https://academic.oup.com/nar/article/43/W1/W39/2467905); [Rawlings et al. 2018](https://academic.oup.com/nar/article/46/D1/D624/4626772)). These patterns, also known as regular expressions (RegExp), were built for each motif:

Aspartic acid (**D**) motif: [VIA][^T][VIFGL][LIVAF]**[D]**[TSADG][GDPS]

Histidine (**H**) motif: **[H]**[GIA][TSDNCM][^R][VCTLIA][AISTG][GSHAL]

Serine (**S**) motif: [G][TN]**[S]**[^A][ASG][STAVCLG][PAG]

Searches in the metagenome based on pattern matching of the regular expression, were performed using the Linux command “grep” typing the following command line in a terminal window:

<pre><code>
$grep -E -B 1 –no-group-separator '[VIA][^T][VIFGL][LIVAF][D][TSADG][GDPS][A-Z]{10,100}[H][GIA][TSDNCM][^R][VCTLIA][AISTG][GSHAL][A-Z]{100,250}[G][TN][S][^A][ASG][STAVCLG][PAG]' lineal_metagenome_aa_file.faa > contigs_with_S8A_proteases.faa
</code></pre>

The previous command line recovers any translated contig that contains the complete sequence of the CT of proteases S8A. Since the contigs may contain more than one single coding sequence (CDS), the following perl command was used to trim the sequences and extract only the S8A conserved region of the CT.

<pre><code>
$perl -lane '{$c=0} if (/>/) {$h=$_; $h=~s/>//} else {while (/([VIA][^T][VIFGL][LIVAF][D][TSADG][GDPS][A-Z]{10,100}[H][GIA][TSDNCM][^R][VCTLIA][AISTG][GSHAL][A-Z]{100,250}[G][TN][S][^A][ASG][STAVCLG][PAG])/g) {{$c++} print ">$c\_$h\n$1"}}' < contigs_with_S8A_proteases.faa > trimmed_S8A_triads.faa
</code></pre>

A bash file including all the commands in sequential order of execution, from translation of contigs to retrieval of the S8A CT subsequences, is provided.
A thorugh explanation of all options used is included (Data_mining_single.sh)

Our methodology is easily scalable; using “for in” loops the commands can be iterated across all files present in a folder effectively mining multiple metagenomes consecutively (Dta_mining_multi.sh)


