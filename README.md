# Modified BLAST for probabilistic sequences
### A local alignment tool for probabilistic DNA sequences

Final Project for Comp 561, Computational Biology Methods and Research (Fall 2019)  
  
The Basic Local Alignment Search Tool (BLAST) is commonly used in bioinformatics to efficiently align and compare a query DNA sequence with a database of sequences. It is available for use online at the [National Center for Biotechnology Information (NCBI) website](https://blast.ncbi.nlm.nih.gov/Blast.cgi), as well as many other sites. The core algorithm is to	find alignments between a short query sequence and a long genome.  
  
However, in certain situations, we may not know exactly the genome of a species, but instead we could have some uncertainty about the identity of nucleotides at a given position (e.g.	at	polymorphic	sites). In that case, the genome of a species may be described using a probability matrix with four rows (A, C, G, T) and L columns, where L is the length of the genome, and the entries are the probabilities of each nucleotide at each position. The goal of this project is to develop an analog of BLAST that works to align a short (non-probabilistic) query sequence against a long probabilistic genome.  
  
While alignment algorithms often uses dynamic programming approach such as the Smith-Waterman algorithm, BLAST combines a heuristic approach with it to largely reduce the complexity of the program, and we modified it for probabilistic DNA sequences with some further optimization. For more details, check out our [report](/BLAST_prob_final_report.pdf) and the [step-by-step demostration](/demo_stepbystep.ipynb). The [demo](/demo.py) shows how to call the pipelined BLAST function with generated query.
