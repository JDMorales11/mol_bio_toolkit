"""ORF (Open Reading Frame) detection in DNA sequences."""

from Bio.Seq import Seq

_STOP_CODONS = {"TAA", "TAG", "TGA"}


def _scan_frame(seq_str: str, frame: int, min_length: int) -> list[dict]:
    """
    Scans a single reading frame for ORFs (ATG...stop codon).

    Returns ORFs with LOCAL coordinates relative to seq_str as given
    (i.e. not yet mapped back to the original input sequence if
    seq_str is a reverse complement). Coordinate mapping is the
    caller's responsibility -- see find_orfs.
    """
    orfs = []
    i = frame
    while i < len(seq_str) - 2:
        codon = seq_str[i:i + 3]
        if codon == "ATG":
            for j in range(i + 3, len(seq_str) - 2, 3):
                stop = seq_str[j:j + 3]
                if stop in _STOP_CODONS:
                    orf_length = j - i
                    if orf_length >= min_length:
                        orf_seq = Seq(seq_str[i:j + 3])
                        orfs.append({
                            "local_start": i,
                            "local_end": j + 3,
                            "length_bp": orf_length,
                            "aa_count": orf_length // 3,
                            "sequence": str(orf_seq),
                            "protein": str(orf_seq[:-3].translate()),
                        })
                    break
        i += 3
    return orfs


def find_orfs(sequence: str | Seq, min_length: int = 100) -> list[dict]:
    """
    Finds all ORFs across all 6 reading frames (3 forward, 3 reverse).

    Genes are not restricted to the strand a sequence happens to be
    written on -- in bacterial genomes and plasmids in particular,
    genes routinely sit on both strands, often overlapping or facing
    each other. Scanning only the forward frames (as in v0.1.x-v0.2.x)
    silently misses any gene encoded on the complementary strand.
    Closes Issue #1.

    Frame numbering follows the standard six-frame translation
    convention: 1, 2, 3 for the forward (+) strand, -1, -2, -3 for
    the reverse (-) strand.

    Coordinate convention: `start` and `end` are always reported
    relative to the ORIGINAL input sequence as given (position 0 =
    first base, 5'->3' as written), regardless of which strand the
    ORF is on. For reverse-strand ORFs this means the raw coordinates
    found on the reverse complement are mapped back via
    `original_start = len(sequence) - local_end` and
    `original_end = len(sequence) - local_start`, so results can be
    sliced directly against the original sequence without extra
    bookkeeping by the caller. The `sequence` and `protein` fields,
    however, always show the ORF as it reads 5'->3' on its own
    strand (i.e. the actual coding sequence that would be
    transcribed/translated), not a slice of the original strand.

    Parameters
    ----------
    sequence : Bio.Seq or str
        Input DNA sequence.
    min_length : int, default 100
    Minimum ORF length in bp, measured from ATG to the base
    immediately before the stop codon (the stop codon itself is
    excluded from the length count but is still included in the
    `sequence` field).

    Returns
    -------
    list of dict
        Each dict has: frame (1..3 or -1..-3), strand ("+"/"-"),
        start, end (original-sequence coordinates), length_bp,
        aa_count, sequence, protein.
    """
    seq_str = str(sequence).upper()
    seq_len = len(seq_str)
    results = []

    # Forward strand: frames 1, 2, 3
    for local_frame in range(3):
        for orf in _scan_frame(seq_str, local_frame, min_length):
            results.append({
                "frame": local_frame + 1,
                "strand": "+",
                "start": orf["local_start"],
                "end": orf["local_end"],
                "length_bp": orf["length_bp"],
                "aa_count": orf["aa_count"],
                "sequence": orf["sequence"],
                "protein": orf["protein"],
            })

    # Reverse strand: frames -1, -2, -3
    rc_str = str(Seq(seq_str).reverse_complement())
    for local_frame in range(3):
        for orf in _scan_frame(rc_str, local_frame, min_length):
            # Map reverse-complement-local coordinates back to the
            # original sequence's coordinate system.
            original_start = seq_len - orf["local_end"]
            original_end = seq_len - orf["local_start"]
            results.append({
                "frame": -(local_frame + 1),
                "strand": "-",
                "start": original_start,
                "end": original_end,
                "length_bp": orf["length_bp"],
                "aa_count": orf["aa_count"],
                "sequence": orf["sequence"],
                "protein": orf["protein"],
            })

    return results