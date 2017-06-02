# MSA_trimmer
Straightforward &amp; minimalistic removal of poorly aligned regions in sequence alignments.

    usage: alignment_trimmer.py [-h] [-p PEP_ALIGNMENT_FASTA]
                                [-c CDS_ALIGNMENT_FASTA] [--trim_gappy TRIM_GAPPY]
                                [--trim_positions TRIM_POSITIONS [TRIM_POSITIONS ...]]
                                [--mask [MASK]]

    A program to remove (trim) or mask columns of fasta multiple sequence
    alignments. Columns are removed if they exceed a user-specified "gappyness"
    treshold, and/or if they are specified by column index. Also works for codon
    alignments. The trimmed alignment is written to [your_alignment]_trimmed.fa

    optional arguments:
    -h, --help            show this help message and exit
    -p PEP_ALIGNMENT_FASTA, --pep_alignment_fasta PEP_ALIGNMENT_FASTA
    -c CDS_ALIGNMENT_FASTA, --cds_alignment_fasta CDS_ALIGNMENT_FASTA
    --trim_gappy TRIM_GAPPY
                            the maximum allowed gap proportion for an MSA column,
                            all columns with a higher proportion of gaps will be
                            trimmed (default: off)
    --trim_positions TRIM_POSITIONS [TRIM_POSITIONS ...]
                            specify columns to be trimmed from the alignment, e.g.
                            "1-10 42" to trim the first ten columns and the 42th
                            column
    --mask [MASK]         mask the alignment instead of trimming. if a character
                            is specified after --mask, that character will be used
                            for masking (default char: "-")