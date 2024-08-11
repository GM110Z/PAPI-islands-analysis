**Heatmap-generator** = Uses PADLOC or AMRFinder outputs to generate heatmaps of presence/absence that can be then used in R

**jarvis.py** = Starting from a list of IDs that where already characterised (by  PADLOC or AMRFinder, or whatever else) jarvis elminates them from the file of representive sequences that is produced by MMseqs to be used for further analysis 

**SPIDERMAN.py** = SPIDERMAN (System for Pathway Identification, Data Extraction, and Retrieval with Marvelous ANalysis) uses a tabled outputted from BlastKoala to map KEGG Ids to Kegg pathways and plots the pathways in a piechart

**edison.py** = Edison (hiddEn moDel proteIn familieS identificatiON)Runs hmmscan vs a PFAM database and processes data to only report hits with e-value <0.01 (and writes them in a tabular format that is better looking than normal pfam)

**calculate-size.sh** =  runs as calculate-size.sh <list-of-nuccore-IDS> Uses efetch from entrez and faidx to calculate size in bp

**plasmid-size-groups.py** = runs as python plasmid-size-groups.py <file_with_sizes>. File with sizes is output from **calculate-size.sh** 
