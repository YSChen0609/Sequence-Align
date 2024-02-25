# Sequence-Align
This is a Repo consists of two modules:

1. A data parser that streams the UNIPROT protein sequence.
2. A Sequence Alignment tool implementing the Needleman-Wunsch algorithm.

## SeqAlign
See `SeqAlign.py`.

Pairwise sequence global alignment to find out their optimal alignment score, optimal alignment, and the corresponding sequence identity.
- Adjustable substitution score matrix (default BLOSUM62)

## CompareProteinSeq
See `CompareProteinSeq.py`.

Fetch and parse the UniProt protein sequence (FASTA format) with this url format: https://rest.uniprot.org/uniprotkb/${seq_id}.fasta.

Experiement the sequence global alignment by calling the `SeqAlign` given a compare set containing the protein sequence name (id).

## Experiments
To get the comparing results, execute `$python __main__.py`.

## Future Work/ Improvements
- Implement pairwise sequence local alignment using Smith–Waterman algorithm.

