###
###   MULTIPLE METAGENOME FILES PROCESSING.
###   BASH file for executing all required steps to recover contigs containing proteins of interest
###   from a ONE or MULTIPLE FASTA formated metagenomic files stored in a single directory.
###

# For this bash file to work in a single run, the files name should be exactly the same.
# The initial file or files should be already assembled in contigs, in FASTA format, and all file names have the ".fna" extension.
#
# Before running check that this file has execcuting permissions.
# Changing permission of this file: $chmod +x Supplementary_file_2.txt
# Running this file: $./Supplementary_file_2.txt


####   1.TRANSLATION   ###
#TRANSLATE the assembled metagenomic contigs into all six posible reading frames


for i in *.fna; do transeq -sequence $i -outseq "aa_"$i".faa" -frame 6 -table  11 -nomethionine; done


# To use "transeq" the EMBOSS package should be installed.

# Parameters and arguments descriptions:
# The expression "for i in *.fna" loops the "transeq" command and applies it sequentially to all files on the working directory with an ".fna" extension.
# Each individual file name is sequentially stored on the "$i" variable.
# "-sequence" indicates the name of the fasta file containing contigs sequences and is automatically retrieved from the "$i" variable.
# "-outseq" indicates name of the new fasta file containing the translations into amino acids of the contigs. It is automatically named prepending the "aa_" tag and the ".faa" extension to the original infile name ("aa_"$i".faa").
# "-frame" indicates the reading frame to be translated. In this case "6" stands for all forward and reverse reading frames.
# "-table" indicates the genetic code (codon table) to be used for translations. In this case "11" stands for "Bacterial" codon table.
# "-nomethionine" indicates that the first codon of the sequence should NOT be translated to methionine and instead translated to its corresponding amino acid. By default transeq treats the first codon, regardless of its sequence, as an start codon and translates it to METHIONINE (M).
# "done" ends the loop when the last file is processed.



###  2.REMOVE LINE BREAKS IN THE TRANSLATED FASTA FILE   ###
# The fasta format  (obtained in the previos step) contains a single-line description followed by lines of sequences data that contains 60 characters in length. 
# The line breaks do not allow the regular expressions to find matching patterns that are divided along multiple lines.
# The following pipeline, that includes "cat" and "perl" commands, removes all line-breaks leaving the full sequence in a single continous line preceded by the single-line description.


for i in aa_*; do cat $i | perl -ne 'if(/^>/){print "\n",$_;next;}else{chomp;print;}' > "lineal_"$i; done



# The "cat" and "perl" commands are installed by default on most linux systems.
# The expression "for i in aa_*" loops the "cat" and "perl" commands and applies them sequentially to all files on the working directory starting with the "aa_" tag.
# The command cat will read the contents of the current file (stored on the "$i" variable) and pass it line by line to the perl command.
# The "e" switch  instructs perl to not look for a script filename on the arguments and instead operate as a one-liner command.
# The "n" switch causes perl to assume a loop around the one-liner command, sequentially applying it to each line feed to it by cat and stores that text in the "$_" variable.
# "if(/^>/)" searches for lines starting with the ">" sign to identify those lines as FASTA headers and adds a new-line (print "\n") at its ends.
# If a line does not start with a ">" sign, it is considered sequence and any newline at its end is deleted. 
# The results will be printed into a new file which name is prepended with the "lineal_" tag.
# "done" ends the loop when the last file is processed.


###   3.REGULAR EXPRESSION MATCHING   ###
###   SEARCH OF TRANSLATED METAGENOME TO FIND CONTIGS CONTAINING THE MOTIF(S) CODED ON A REGULAR EXPRESSION   ###


for i in lineal_*; do grep -E -B 1 --no-group-separator '[VIA][^T][VIFGL][LIVAF][D][TSADG][GDPS][A-Z]{10,100}[H][GIA][TSDNCM][^R][VCTLIA][AISTG][GSHAL][A-Z]{100,250}[G][TN][S][^A][ASG][STAVCLG][PAG]' $i > "S8_contigs_"$i; done


# The expression "for i in lineal_*" loops the "grep" command and applies it sequentially to all files on the working directory starting with the "lineal_" tag.
# The "grep" command is installed by default on most linux systems.
# "-E" flag forces grep to use the extended regular expression syntaxis so "{" and "}" are considered metacharacters and not matched literally.
# "-B 1" prints one (1) preceding line (the fasta header) in addition to the line containing the matched regular expression (in this case a contig containing a full S8 triad).
# "--no-group-separator" indicates to not print a new line with two dashes (--). By default grep will separate each match found from each other with a new line with two dashes (--).
# The regular expression is surrounded by single quotes ('). 
# "lineal_metagenome_aa_file.faa" is the input fasta file with the translated and linearized metagenomic contigs. 
# The results will be printed into a new file which name is prepended with the "S8_contigs_" tag. 
# "done" ends the loop when the last file is processed.


###   4.TRIMMING OF THE CONTIGS CONTAINING MATCHING FULL REGULAR EXPRESSION   ### 
# The following perl command is to trim the contig sequences and extract all subsequences correspoding to the Catalitic Triad.


for i in S8_contigs_*; do perl -lane '{$c = 0} if (/>/) {$h=$_; $h=~s/>//} else { while (/([VIA][^T][VIFGL][LIVAF][D][TSADG][GDPS][A-Z]{10,100}[H][GIA][TSDNCM][^R][VCTLIA][AISTG][GSHAL][A-Z]{100,250}[G][TN][S][^A][ASG][STAVCLG][PAG])/g) { {$c++} print ">$c\_$h\n$1"}}' < $i > "S8_triad_"$i; done


# The expression "for i in S8_contigs_*" loops the "perl" command and applies it sequentially to all files on the working directory starting with the "S8_contigs_" tag.
# The "e" switch  instructs perl to not look for a script filename on the arguments and instead operate as a one-liner command.
# The "n" switch causes perl to assume a loop around the one-liner command.
# The "l" switch enables automatic line-ending processing.
# The "a" switch turns on autosplit mode when used with a "n".
# "{$c = 0}" initialize a counter to add to each FASTA header to avoid repetitions when multiple CT are found inside the same contig."{$c++}" increases the counter one unit after each positive match of a CT in a contig.
# "(/>/) {$h=$_; $h=~s/>//}" identifies FASTA headers (lines containing a ">" sign) and rewrites those lines adding the counter.
# If a ">" is not found, the regular expression is searched on the line. Whenever it is matched, it is printed in a new-line after its corresponding FASTA header. The "g" switch after the regular expressions matches it "globally" continuing to search for it even if it has a match on the contig.
# The input fasta file name is automatically retrieved from the "$i" variable. 
# The results will be printed into a new file which name is prepended with the "S8_triad_" tag. 
# "done" ends the loop when the last file is processed.

