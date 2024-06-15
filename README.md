# capture_mock_communities
The python and r-markdown script from the paper: "Enriching barcoding markers in environmental samples utilising a phylogenetic probe design: insights from mock communities"

The raw data required for running these scripts is accessible from Figshare. The Python script 'mapped_read_parser_V1.py' can be used on any Bowtie2 BAM file where reads are mapped locally to bait sequences. It calculates the GC-content of the longest ungapped sequence and the number of mismatches in this part of the alignment.

The 'chimeric_seq_finder.py' script examines local alignments to identify reads that are mapped to two reference sequences, with one forward match and one reverse match, which would indicate a chimeric sequence. 
