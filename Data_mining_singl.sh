###
###   SINGLE METAGENOME FILE PROCESSING.
###   BASH file for executing all required steps to recover contigs containing proteins of interest
###   from a single FASTA formated metagenomic file.
###

# For this bash file to work in a single run, the files name should be exactly the same.
# The initial file should be already assembled in contigs, in FASTA format, and named "metagenome_contigs.fna".
#
# Before running check that this file has execcuting permissions.
# Changing permission of this file: $chmod +x Supplementary_file_1.txt
# Running this file: $./Supplementary_file_1.txt


####   1.TRANSLATION   ###
#TRANSLATE the assembled metagenomic contigs into all six posible reading frames


transeq -sequence metagenome_contigs.fna -outseq metagenome_aa_file.faa -frame 6 -table  11 -nomethionine


# To use "transeq" the EMBOSS package should be installed.
# Parameters and arguments descriptions:
# "-sequence" indicates the fasta file containing contigs sequences. Change the name of the fasta file accordingly.
# "-outseq" indicates name (and path) of the new fasta file containing the translations into amino acids of the contigs. Change name of output file accordingly.
# "-frame" indicates the reading frame to be translated. In this case "6" stands for all forward and reverse reading frames.
# "-table" indicates the genetic code (codon table) to be used for translations. In this case "11" stands for "Bacterial" codon table.
# "-nomethionine" indicates that the first codon of the sequence should NOT be translated to methionine and instead translated to its corresponding amino acid. By default transeq treats the first codon, regardless of its sequence, as an start codon and translates it to METHIONINE (M).



###  2.REMOVE LINE BREAKS IN THE TRANSLATED FASTA FILE   ###
# The fasta format  (obtained in the previos step) contains a single-line description followed by lines of sequences data that contains 60 characters in length. 
# The line breaks do not allow the regular expressions to find matching patterns that are divided along multiple lines.
# The following pipeline, that includes "cat" and "perl" commands, removes all line-breaks leaving the full sequence in a single continous line preceded by the single-line description.


cat metagenome_aa_file.faa | perl -ne 'if(/^>/){print "\n",$_;next;}else{chomp;print;}' > lineal_metagenome_aa_file.faa



# The "cat" and "perl" commands are installed by default on most linux systems.
# The command cat will read the contents of the metagenome_aa_file.faa file and pass it line by line to the perl command.
# The "e" switch  instructs perl to not look for a script filename on the arguments and instead operate as a one-liner command.
# The "n" switch causes perl to assume a loop around the one-liner command, sequentially applying it to each line feed to it by cat and stores that text in the "$_" variable.
# "if(/^>/)" searches for lines starting with the ">" sign to identify those lines as FASTA headers and adds a new-line (print "\n") at its ends.
# If a line does not start with a ">" sign, it is considered sequence and any newline at its end is deleted. 
# The results will be printed into a new file named "lineal_metagenome_aa_file.faa".



###   3.REGULAR EXPRESSION MATCHING   ###
###   SEARCH OF TRANSLATED METAGENOME TO FIND CONTIGS CONTAINING THE MOTIF(S) CODED ON A REGULAR EXPRESSION   ###


grep -E -B 1 --no-group-separator '[VIA][^T][VIFGL][LIVAF][D][TSADG][GDPS][A-Z]{10,100}[H][GIA][TSDNCM][^R][VCTLIA][AISTG][GSHAL][A-Z]{100,250}[G][TN][S][^A][ASG][STAVCLG][PAG]' lineal_metagenome_aa_file.faa > contigs_with_S8_proteases.faa


# The "grep" command is installed by default on most linux systems.
# "-E" flag forces grep to use the extended regular expression syntaxis so "{" and "}" are considered metacharacters and not matched literally.
# "-B 1" prints one (1) preceding line (the fasta header) in addition to the line containing the matched regular expression (in this case a contig containing a full S8 triad).
# "--no-group-separator" indicates to not print a new line with two dashes (--). By default grep will separate each match found from each other with a new line with two dashes (--).
# The regular expression is surrounded by single quotes ('). 
# "lineal_metagenome_aa_file.faa" is the input fasta file with the translated and linearized metagenomic contigs. 
# The results will be printed into a new file named "contigs_with_S8_proteases.faa". 



###   4.TRIMMING OF THE CONTIGS CONTAINING MATCHING FULL REGULAR EXPRESSION   ### 
# The following perl command is to trim the contig sequences and extract all subsequences correspoding to the Catalitic Triad.


perl -lane '{$c = 0} if (/>/) {$h=$_; $h=~s/>//} else { while (/([VIA][^T][VIFGL][LIVAF][D][TSADG][GDPS][A-Z]{10,100}[H][GIA][TSDNCM][^R][VCTLIA][AISTG][GSHAL][A-Z]{100,250}[G][TN][S][^A][ASG][STAVCLG][PAG])/g) { {$c++} print ">$c\_$h\n$1"}}' < contigs_with_S8_proteases.faa > S8_triad.faa


# The "e" switch  instructs perl to not look for a script filename on the arguments and instead operate as a one-liner command.
# The "n" switch causes perl to assume a loop around the one-liner command.
# The "l" switch enables automatic line-ending processing.
# The "a" switch turns on autosplit mode when used with a "n".
# "{$c = 0}" initialize a counter to add to each FASTA header to avoid repetitions when multiple CT are found inside the same contig."{$c++}" increases the counter one unit after each positive match of a CT in a contig.
# "(/>/) {$h=$_; $h=~s/>//}" identifies FASTA headers (lines containing a ">" sign) and rewrites those lines adding the counter.
# If a ">" is not found, the regular expression is searched on the line. Whenever it is matched, it is printed in a new-line after its corresponding FASTA header. The "g" switch after the regular expressions matches it "globally" continuing to search for it even if it has a match on the contig.
# "contigs_with_S8_proteases.faa" is the input fasta file. 
# The results will be printed into a new file named "S8_triad.faa"



