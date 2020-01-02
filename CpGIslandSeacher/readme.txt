README (added at UCSC, update by Jian-Hong Sun)

Step 1 Compile:
"make" --> gcc readseq.c cpg_lh_new.c -o cpglh_new.exe

Step 2 run:
for sigle inputfile:
	run: >cpgTest.exe inputfile(*.fasta format)
	Note: The output directory of this program has been set to "D: \ resutls \"
	         (you can change the 109 lines of cpg_lh_new.c)

for batch file processing(highly recommended):
	use .bat file (E.g: CpGIslandSeacher.bat , E.g: CpGIslandSeacher2.bat, E.g: CpGIslandSeacher3.bat)
	Note: The list of file names written in "name.txt", all the inputfile are placed in same folder with the program

############################################################################################
We've been running this on hard-masked sequence (RepeatMasker and TRF 
with period <= 12 results are masked to 'N') in order to avoid CpG 
islands in Alu repeats in human.

The original cpg.c was written by Gos Miklem from the Sanger Center.  
LaDeana Hillier added some modifications --> cpg_lh.c, and UCSC has 
added some further modifications to cpg_lh.c, and Jian-Hong Sun 
has added some further modifications  --> cpg_lh_new.c, so that its expected 
number of CpGs in an island is calculated as described in 
  Gardiner-Garden, M. and M. Frommer, 1987 
  CpG islands in vertebrate genomes. J. Mol. Biol. 196:261-282.

    Expected = (Number of C's * Number of G's) / Length

Instead of a sliding-window search for CpG islands, this cpg program
uses a running-sum score where a 'C' followed by a 'G' increases the
score by 17 and anything else decreases the score by 1.  When the
score transitions from positive to 0 (and at the end of the sequence),
the sequence in the current span is evaluated to see if it qualifies
as a CpG island (>200 bp length, >50% GC, >0.6 ratio of observed CpG
to expected).  Then the search recurses on the span from the position
with the max running score up to the current position.
