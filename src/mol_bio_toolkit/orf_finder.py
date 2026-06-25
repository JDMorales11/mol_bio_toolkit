"""ORF (Open Reading Frame) detection in DNA sequences."""

from Bio.Seq import Seq

def find_orfs(sequence, min_length=100):
    """
    Encuentra todos los ORFs en los 3 marcos de lectura forward.

    Parameters
    ----------
    sequence : Bio.Seq o str - secuencia de ADN
    min_length : int - longitud mínima del ORF en bp (sin contar el STOP codon)

    Returns
    -------
    list of dict: frame, start, end, lenght_bp, aa_count, sequence, protein
    """
    orfs = []
    seq_str = str(sequence).upper()

    for frame in range(3):
        i = frame
        while i < len(seq_str) - 2:
            codon = seq_str[i:i+3]
            if codon == "ATG": 
                for j in range(i+3, len(seq_str)-2, 3):
                    stop = seq_str[j:j+3]
                    if stop in ["TAA", "TAG", "TGA"]:
                        orf_length = j - i
                        if orf_length >= min_length:
                            orf_seq = Seq(seq_str[i:j+3])
                            orfs.append({
                                "frame": frame + 1,
                                "start": i,
                                "end": j + 3,
                                "length_bp": orf_length,
                                "aa_count": orf_length // 3,
                                "sequence": str(orf_seq),
                                "protein": str(orf_seq[:-3].translate()),
                            })
                        break
            i += 3
    return orfs
