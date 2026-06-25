"""Optimizacion de codones para expresion heterologa."""

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

# Tabla de codones preferidos por organismo.
# Fuente: Codon Usage Database (Nakamura et al. 2000).
# Cada valor es el codon sinonimo de mayor frecuencia relativa
# para ese aminoacido en el genoma de referencia del organismo.
CODON_TABLES = {
    "e_coli_k12": {
        "A": "GCT", "R": "CGT", "N": "AAC", "D": "GAT", "C": "TGC",
        "Q": "CAG", "E": "GAA", "G": "GGT", "H": "CAC", "I": "ATT",
        "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCG",
        "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAT", "V": "GTT",
        "*": "TAA",
    },
    "s_cerevisiae": {
        "A": "GCT", "R": "AGA", "N": "AAC", "D": "GAT", "C": "TGT",
        "Q": "CAA", "E": "GAA", "G": "GGT", "H": "CAC", "I": "ATT",
        "L": "TTG", "K": "AAA", "M": "ATG", "F": "TTC", "P": "CCA",
        "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAC", "V": "GTT",
        "*": "TAA",
    },
}


def codon_optimize(cds_sequence, host="e_coli_k12"):
    """
    Recodifica un CDS usando los codones preferidos del organismo host.

    El codigo genetico es degenerado: 61 codones codifican solo 20
    aminoacidos, asi que la mayoria de aminoacidos tienen multiples
    codones sinonimos posibles. Cada organismo usa esos sinonimos con
    distinta frecuencia (sesgo de uso de codones), reflejando la
    abundancia relativa de tRNAs en esa celula. Un gen mal optimizado
    para su huesped puede traducirse lento o de forma ineficiente.

    Esta funcion preserva la proteina resultante de forma exacta:
    solo cambia los codones, nunca la secuencia de aminoacidos.

    Parameters
    ----------
    cds_sequence : Bio.Seq or str
        CDS completo en marco de lectura (divisible por 3, incluye
        codon de stop).
    host : str, default "e_coli_k12"
        Organismo objetivo. Opciones disponibles en CODON_TABLES:
        "e_coli_k12", "s_cerevisiae".

    Returns
    -------
    dict
        Incluye: host, original_seq, optimized_seq, aa_length,
        codon_changes, identity_%, original_gc%, optimized_gc%,
        y codon_log (DataFrame con detalle posicion a posicion).

    Raises
    ------
    KeyError
        Si host no esta en CODON_TABLES.
    """
    if host not in CODON_TABLES:
        raise KeyError(
            f"Host '{host}' no reconocido. Opciones: {list(CODON_TABLES)}"
        )
    codon_table = CODON_TABLES[host]

    seq_str = str(cds_sequence).upper()
    protein = str(Seq(seq_str[:-3]).translate())

    optimized = ""
    changes = 0
    codon_log = []

    for i, aa in enumerate(protein):
        original_codon = seq_str[i * 3: i * 3 + 3]
        preferred_codon = codon_table.get(aa, original_codon)
        optimized += preferred_codon
        changed = original_codon != preferred_codon
        if changed:
            changes += 1
        codon_log.append({
            "position": i + 1,
            "aa": aa,
            "original": original_codon,
            "optimized": preferred_codon,
            "changed": changed,
        })

    optimized += codon_table.get("*", "TAA")
    identity = (len(protein) + 1 - changes) / (len(protein) + 1) * 100

    return {
        "host": host,
        "original_seq": seq_str,
        "optimized_seq": optimized,
        "aa_length": len(protein),
        "codon_changes": changes,
        "identity_%": round(identity, 1),
        "original_gc%": round(gc_fraction(Seq(seq_str)) * 100, 1),
        "optimized_gc%": round(gc_fraction(Seq(optimized)) * 100, 1),
        "codon_log": pd.DataFrame(codon_log),
    }