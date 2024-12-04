**Heatmap-generator** = Uses PADLOC or AMRFinder outputs to generate heatmaps of presence/absence that can be then used in R

**jarvis.py** = Filters blast ouput to produce a table of regions with start and stop that can be fed to efetch. It solves the problem of local alignments by separating regions that are more than a certain distance (i.e. 180k for PAPi-1). It then calculate and groups regions by size intervals

**SPIDERMAN.py** = SPIDERMAN (System for Pathway Identification, Data Extraction, and Retrieval with Marvelous ANalysis) uses a tabled outputted from BlastKoala to map KEGG Ids to Kegg pathways and plots the pathways in a piechart

**edison.py** = Edison (hiddEn moDel proteIn familieS identificatiON)Runs hmmscan vs a PFAM database and processes data to only report hits with e-value <0.01 (and writes them in a tabular format that is better looking than normal pfam)

**Prokkaloop** : Edits headers of fasta files when downloading slices and then runs prokka



Order to run scripts

1)Jarvis
2)Download with parallel efetch and then run prokkaloop
