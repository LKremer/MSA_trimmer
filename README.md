# MSA_trimmer: Straightforward &amp; minimalistic removal of poorly aligned regions in **M**ultiple **S**equence **A**lignments
This is a program to remove (trim) or mask columns of multiple sequence alignments in fasta format.
Columns are removed if they exceed a user-specified "gappyness" treshold, i.e. a threshold of
0.5 (``--trim_gappy 0.5``) would trim all columns where more than half of the sequences have a gap. You may also "hand-pick" columns
or whole regions of the alignment to be trimmed (``--trim_positions``).

If you want to preserve the original coordinates of the alignment, you may also chose to "mask" your alignment
instead of trimming it. This will not remove any columns, instead the whole column will be replaced with
a character of your choice (the gap character "-" by default).

Trimming works for both nucleotide/codon/CDS (``-c``) and peptide (``-p``) alignments. You may even trim a pair
of corresponding CDS and peptide alignments at the same time, to ensure the trimmed alignments will match. In that
case, only whole codons (chunks of three bases each) will be removed in order not to introduce a frame shift.
This is useful for applications that need a codon alignment such as CODEML.

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
