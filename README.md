**Heatmap-generator** = Uses PADLOC or AMRFinder outputs to generate heatmaps of presence/absence that can be then used in R

**jarvis.py** = Filters blast ouput to produce a table of regions with start and stop that can be fed to efetch. It solves the problem of local alignments by separating regions that are more than a certain distance (i.e. 180k for PAPi-1). It then calculate and groups regions by size intervals

**SPIDERMAN.py** = SPIDERMAN (System for Pathway Identification, Data Extraction, and Retrieval with Marvelous ANalysis) uses a tabled outputted from BlastKoala to map KEGG Ids to Kegg pathways and plots the pathways in a piechart

**edison.py** = Edison (hiddEn moDel proteIn familieS identificatiON)Runs hmmscan vs a PFAM database and processes data to only report hits with e-value <0.01 (and writes them in a tabular format that is better looking than normal pfam)

**Prokkaloop** : Edits headers of fasta files when downloading slices and then runs prokka

**split-Mmseqs-representative-fastas.sh**: After running mmseq to deruplicate, it splits the representative sequence files from multifastas to single fasta that are compatible with prokkaloop script

**PHORIFIC.py**: Needs a tab separated file from MMseqs clustering, and genbank files of interest to be run.  

**SantasHelper.py** Combines PHORIFIC output with PADLOC, AMRFInder and Antidefensefinder

**Clip.py** Combines Probe and Jarvis outputs to ensure that probe only picks regions included in Jarvis output

**Probe.py** uses hmm models to find co-localised genes. 

**fasta-header-extarctor.sh** From a multifasta file produced by mmseqs, it extracts nuccore IDs and coordinates from headers of representative seq file for further analysis

**representative.hotspot.selection.py**: To use with bigger dataset. Starts from a list of nuccore ID you can get from padloc or defense finder. The it compares it with the whole blast output of Jarvis to select representative genomes that have defense systems and can be used to analyse to locate hotspot boundaries

**heatmap-plots.py**: For those islands with multiple hotspots, produces a presence/absence matrix of each hotspot for each nuccore ID (clustered as a dendrogram based on hotspot presence/absence) 

**detective.py 	DETECT Ice oVErlap** : It compares files with jarvis coordinates for each pathogenicity islands and detects for overlap. It produces a table for coordinates and size of overlap and a scatter plot to summarise this

**viper.py  – Virulence Identification and Protein Elimination via Refinement**: It needs a protein file.  WIll run PSI-BLAST using the inputed proteome vs a database of your choice (e.g. VFDB). 
*You have to make your own makeblastdb database and then can input the name as sys.argv[1]*
*proteome_file = sys.argv[2]*
*output_dir = sys.argv[3]*


