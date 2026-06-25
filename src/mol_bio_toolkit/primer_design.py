"""Diseno de primers para PCR: Tm, GC-clamp, deteccion de hairpins."""

from Bio.Seq import Seq


def calculate_tm(primer):
    """
    Calcula la temperatura de melting (Tm) de un primer.

    Usa la regla de Wallace para oligos cortos (< 14 nt), donde cada
    par de bases contribuye un valor fijo de energia de hibridacion.
    Para primers mas largos (>= 14 nt) usa una formula con correccion
    por longitud total, porque el efecto proporcional del GC% domina
    sobre el conteo absoluto de bases.

    Parameters
    ----------
    primer : str
        Secuencia del primer (5' -> 3').

    Returns
    -------
    float
        Tm estimada en grados Celsius.
    """
    s = primer.upper()
    gc = s.count("G") + s.count("C")
    at = s.count("A") + s.count("T")
    n = len(s)
    if n < 14:
        return 4 * gc + 2 * at
    return 64.9 + 41 * (gc - 16.4) / n


def gc_clamp(primer, n=3):
    """
    Verifica si el primer tiene GC-clamp en su extremo 3'.

    La DNA polimerasa extiende el primer desde su extremo 3'. Un G o C
    ahi (3 puentes de hidrogeno vs 2 de A-T) da un punto de union mas
    estable, reduciendo el riesgo de que la polimerasa se desprenda
    antes de iniciar la extension.

    Parameters
    ----------
    primer : str
    n : int, default 3
        Cuantos nucleotidos del extremo 3' revisar.

    Returns
    -------
    bool
    """
    tail = primer[-n:].upper()
    return any(b in "GC" for b in tail)


def has_hairpin(primer, min_stem=4):
    """
    Deteccion heuristica de estructura secundaria tipo hairpin.

    Busca si algun segmento del primer es complementario inverso a
    otra region del mismo primer. Si existe, esas dos regiones pueden
    aparearse entre si en vez de unirse a la secuencia template,
    secuestrando el primer y reduciendo la eficiencia de la PCR.

    Nota: esta es una heuristica simple basada en coincidencia de
    secuencia, no un calculo termodinamico de energia libre como el
    que usan herramientas como Primer3. Sirve como filtro rapido,
    no como sustituto de validacion experimental.

    Parameters
    ----------
    primer : str
    min_stem : int, default 4
        Longitud minima del segmento complementario para considerarlo
        un hairpin potencial.

    Returns
    -------
    bool
    """
    s = primer.upper()
    rc = str(Seq(s).reverse_complement())
    for i in range(len(s) - min_stem * 2):
        stem = s[i:i + min_stem]
        if stem in rc:
            return True
    return False


def _evaluate_primer(seq, tm_min, tm_max):
    """Helper interno: calcula todas las metricas de un primer."""
    tm = calculate_tm(seq)
    gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
    return {
        "sequence": seq,
        "length": len(seq),
        "Tm_C": round(tm, 1),
        "GC_%": round(gc, 1),
        "GC_clamp": gc_clamp(seq),
        "hairpin": has_hairpin(seq),
        "Tm_ok": tm_min <= tm <= tm_max,
    }


def find_optimal_primer(sequence, anchor, direction="forward",
                        len_range=(18, 30), tm_min=55, tm_max=65):
    """
    Busca iterativamente la longitud de primer que cae en el rango de Tm.

    En vez de fijar una longitud arbitraria (ej. 20 nt) que puede o no
    caer en el rango de Tm deseado segun el GC% local de la secuencia,
    esta funcion prueba todas las longitudes posibles en len_range y
    devuelve la mejor opcion encontrada: prioriza primers con Tm dentro
    del rango, GC-clamp presente y sin hairpin. Si ninguna longitud
    cumple las tres condiciones, devuelve la que mas cerca este del
    centro del rango de Tm.

    Parameters
    ----------
    sequence : Bio.Seq or str
        Secuencia completa de la cual extraer el primer.
    anchor : int
        Posicion de anclaje. Para forward, el primer empieza aqui.
        Para reverse, el primer termina aqui (se toma reverse complement).
    direction : {"forward", "reverse"}, default "forward"
    len_range : tuple(int, int), default (18, 30)
        Rango de longitudes de primer a evaluar.
    tm_min, tm_max : float
        Rango de Tm aceptable en grados Celsius.

    Returns
    -------
    dict
        Metricas del mejor primer encontrado, mas la lista completa de
        candidatos evaluados bajo la key "all_candidates".

    Raises
    ------
    ValueError
        Si direction no es "forward" o "reverse", o si no se puede
        generar ningun candidato con los parametros dados.
    """
    seq_str = str(sequence).upper()
    candidates = []

    for length in range(len_range[0], len_range[1] + 1):
        if direction == "forward":
            raw = seq_str[anchor: anchor + length]
        elif direction == "reverse":
            raw = seq_str[anchor - length: anchor]
            raw = str(Seq(raw).reverse_complement())
        else:
            raise ValueError("direction debe ser 'forward' o 'reverse'")

        if len(raw) < length:
            continue  # se salio de los limites de la secuencia

        candidates.append(_evaluate_primer(raw, tm_min, tm_max))

    if not candidates:
        raise ValueError(
            "No se pudo generar ningun primer candidato; revisa "
            "que 'anchor' y 'len_range' esten dentro de la secuencia."
        )

    # Prioridad: Tm_ok + GC_clamp + sin hairpin > solo Tm_ok > mas cercano al centro
    tm_center = (tm_min + tm_max) / 2

    def score(c):
        perfect = c["Tm_ok"] and c["GC_clamp"] and not c["hairpin"]
        return (
            0 if perfect else (1 if c["Tm_ok"] else 2),
            abs(c["Tm_C"] - tm_center),
        )

    best = min(candidates, key=score)
    best["all_candidates"] = candidates
    return best


def design_primers(sequence, amplicon_start, amplicon_end,
                   len_range=(18, 30), tm_min=55, tm_max=65):
    """
    Disena un par forward/reverse optimizado para amplificar una region.

    Cada primer del par se busca de forma independiente dentro de
    len_range para maximizar la probabilidad de caer en el rango de
    Tm deseado, en vez de usar una longitud fija arbitraria.

    Parameters
    ----------
    sequence : Bio.Seq or str
    amplicon_start, amplicon_end : int
        Coordenadas de la region a amplificar (0-indexed).
    len_range : tuple(int, int), default (18, 30)
    tm_min, tm_max : float

    Returns
    -------
    dict
        Keys: "forward", "reverse" (cada uno con sus metricas) y
        "amplicon_size".
    """
    forward = find_optimal_primer(
        sequence, anchor=amplicon_start, direction="forward",
        len_range=len_range, tm_min=tm_min, tm_max=tm_max,
    )
    reverse = find_optimal_primer(
        sequence, anchor=amplicon_end, direction="reverse",
        len_range=len_range, tm_min=tm_min, tm_max=tm_max,
    )
    return {
        "forward": forward,
        "reverse": reverse,
        "amplicon_size": amplicon_end - amplicon_start,
    }