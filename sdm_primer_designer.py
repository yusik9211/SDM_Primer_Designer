#!/usr/bin/env python3
"""
SDM (Site-Directed Mutagenesis) Double Primer Designer

Design rules enforced:
  - Mutation included in both primers
  - Mutation >= 4 nt from 5' terminus
  - Mutation >= 8 nt from 3' terminus
  - >= 8 non-overlapping bases at 3' end of each primer
  - 3' terminal base is G or C
  - Tm >= 78 C  (Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch)

Amino acid mode:
  - Specify protein position + original AA + new AA
  - Best codon selected by E. coli K-12 codon usage frequency
  - Primers named {Gene}_{Mut}_F / {Gene}_{Mut}_R (e.g. MaEgtB_G100A_F)

Copyright (c) 2026 Yusik Kim, Korea University
Contact: yusik9211@korea.ac.kr

Reference:
  Zheng, L., U. Baumann, and Jean-Louis Reymond. (2004).
  An efficient one-step site-directed and site-saturation mutagenesis protocol.
  Nucleic Acids Res. 32(14): e115.
"""

import sys

_COMP = str.maketrans("ATGCatgcNn", "TACGtacgNn")

# --------------------------- codon tables ----------------------------

CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# E. coli K-12 codon usage (per 1000 codons)
# Source: Codon Usage Database (Kazusa), E. coli K-12 MG1655
ECOLI_CODON_FREQ = {
    "TTT": 22.0, "TTC": 17.4,
    "TTA": 14.1, "TTG": 13.6, "CTT": 11.8, "CTC": 10.4, "CTA":  3.6, "CTG": 52.4,
    "ATT": 29.0, "ATC": 25.1, "ATA":  7.5,
    "ATG": 27.4,
    "GTT": 18.1, "GTC": 15.3, "GTA": 10.6, "GTG": 26.1,
    "TCT":  9.0, "TCC":  8.7, "TCA":  7.2, "TCG":  8.5, "AGT":  8.7, "AGC": 16.0,
    "CCT":  7.3, "CCC":  5.0, "CCA":  8.5, "CCG": 23.3,
    "ACT":  9.5, "ACC": 23.3, "ACA":  7.5, "ACG": 14.1,
    "GCT": 15.4, "GCC": 26.1, "GCA": 21.1, "GCG": 33.5,
    "TAT": 16.3, "TAC": 12.3,
    "TAA":  2.0, "TAG":  0.3, "TGA":  1.1,
    "CAT": 13.3, "CAC":  9.0,
    "CAA": 15.4, "CAG": 28.9,
    "AAT": 18.4, "AAC": 21.5,
    "AAA": 33.6, "AAG": 10.2,
    "GAT": 31.5, "GAC": 19.2,
    "GAA": 39.1, "GAG": 17.8,
    "TGT":  5.0, "TGC":  6.3,
    "TGG": 15.3,
    "CGT": 20.6, "CGC": 22.3, "CGA":  3.6, "CGG":  5.4, "AGA":  2.1, "AGG":  1.2,
    "GGT": 24.4, "GGC": 29.2, "GGA":  8.0, "GGG": 11.0,
}

AA_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Stop",
}

# --------------------------- codon helpers ---------------------------

def translate_codon(codon):
    return CODON_TABLE.get(codon.upper(), "?")


def ecoli_codons_ranked(aa):
    """All codons for aa, sorted by E. coli K-12 frequency descending."""
    aa = aa.upper()
    return sorted(
        [(c, f) for c, f in ECOLI_CODON_FREQ.items() if CODON_TABLE.get(c) == aa],
        key=lambda x: x[1],
        reverse=True,
    )


def best_ecoli_codon(aa):
    ranked = ecoli_codons_ranked(aa)
    if not ranked:
        raise ValueError("Unknown amino acid: '%s'" % aa)
    return ranked[0][0]


def aa_mutation_to_dna(template, aa_pos, orig_aa, new_aa, cds_start=1):
    T = validate_dna(template)
    orig_aa = orig_aa.upper()
    new_aa  = new_aa.upper()

    if orig_aa not in AA_NAMES:
        raise ValueError("Unknown original amino acid: '%s'" % orig_aa)
    if new_aa not in AA_NAMES:
        raise ValueError("Unknown new amino acid: '%s'" % new_aa)
    if new_aa == "*":
        raise ValueError("Stop codon as target is not supported in this tool")

    nt_pos = cds_start + (aa_pos - 1) * 3
    cs = nt_pos - 1
    ce = cs + 3

    if ce > len(T):
        raise ValueError(
            "AA position %d (nt %d-%d) exceeds template length %d"
            % (aa_pos, nt_pos, nt_pos + 2, len(T))
        )

    orig_codon = T[cs:ce]
    tmpl_aa = translate_codon(orig_codon)

    if tmpl_aa != orig_aa:
        raise ValueError(
            "AA pos %d: template codon '%s' = %s (%s), expected %s (%s)"
            % (aa_pos, orig_codon, tmpl_aa, AA_NAMES.get(tmpl_aa, "?"),
               orig_aa, AA_NAMES.get(orig_aa, "?"))
        )

    new_codon = best_ecoli_codon(new_aa)

    info = {
        "aa_pos":     aa_pos,
        "orig_aa":    orig_aa,
        "new_aa":     new_aa,
        "orig_codon": orig_codon,
        "new_codon":  new_codon,
        "ranked":     ecoli_codons_ranked(new_aa),
        "nt_pos":     nt_pos,
        "cds_start":  cds_start,
    }
    return nt_pos, orig_codon, new_codon, info


# --------------------------- DNA core --------------------------------

def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def validate_dna(seq):
    seq = seq.upper().replace(" ", "").replace("\n", "").replace("\t", "")
    bad = set(seq) - set("ATGCN")
    if bad:
        raise ValueError("Invalid DNA characters: " + ", ".join(sorted(bad)))
    return seq


def calc_tm(seq, n_mm=0):
    """Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch"""
    N = len(seq)
    if N == 0:
        return float("-inf")
    gc = sum(1 for b in seq.upper() if b in "GC")
    return 81.5 + 0.41 * (gc * 100.0 / N) - 675.0 / N - (n_mm * 100.0 / N)


def gc_pct(seq):
    if not seq:
        return 0.0
    return sum(1 for b in seq.upper() if b in "GC") / len(seq) * 100.0


def count_mm(orig, new):
    n = min(len(orig), len(new))
    subs = sum(1 for a, b in zip(orig.upper(), new.upper()) if a != b)
    return subs + abs(len(orig) - len(new))


def _primer_candidates(tmpl, ms, me, mut, n_mm, fwd, max_len=70):
    ML5, ML3 = 4, 8
    results = []
    max_fl5 = min(ms, max_len)
    max_fl3 = min(len(tmpl) - me if fwd else ms, max_len)

    for fl5 in range(ML5, max_fl5 + 1):
        for fl3 in range(ML3, max_fl3 + 1):
            if fwd:
                s, e = ms - fl5, me + fl3
            else:
                s, e = ms - fl3, me + fl5

            if s < 0 or e > len(tmpl):
                continue

            region = tmpl[s:ms] + mut + tmpl[me:e]
            seq = region if fwd else revcomp(region)

            if seq[-1] not in "GC":
                continue

            tm = calc_tm(seq, n_mm)
            if tm >= 78.0:
                results.append({
                    "seq": seq,
                    "tm":  round(tm, 2),
                    "len": len(seq),
                    "gc":  round(gc_pct(seq), 1),
                    "fl5": fl5,
                    "fl3": fl3,
                    "ts":  s + 1,
                    "te":  e,
                })

    results.sort(key=lambda x: x["len"])
    return results


def find_primers(tmpl, ms, me, mut, n_mm, max_len=70):
    fwd_list = _primer_candidates(tmpl, ms, me, mut, n_mm, fwd=True,  max_len=max_len)
    rev_list = _primer_candidates(tmpl, ms, me, mut, n_mm, fwd=False, max_len=max_len)

    best_fwd = best_rev = None
    best_total = float("inf")

    for f in fwd_list:
        if f["len"] + (rev_list[0]["len"] if rev_list else 0) >= best_total:
            continue
        for r in rev_list:
            if f["len"] + r["len"] >= best_total:
                break
            if f["fl3"] < r["fl5"] + 8:
                continue
            if r["fl3"] < f["fl5"] + 8:
                continue
            best_total = f["len"] + r["len"]
            best_fwd, best_rev = f, r

    return best_fwd, best_rev


def design(template, pos, orig, mut):
    T = validate_dna(template)
    O = validate_dna(orig) if orig else ""
    M = validate_dna(mut)  if mut  else ""

    ms = pos - 1
    me = ms + len(O)

    if me > len(T):
        raise ValueError(
            "Position %d + length %d exceeds template (%d bp)" % (pos, len(O), len(T))
        )
    if O and T[ms:me] != O:
        raise ValueError(
            "Template[%d:%d] = '%s'  but original given = '%s'"
            % (pos, pos + len(O) - 1, T[ms:me], O)
        )

    n_mm = count_mm(O, M)
    fwd, rev = find_primers(T, ms, me, M, n_mm)

    fwd_nonoverlap = (fwd["fl3"] - rev["fl5"]) if (fwd and rev) else 0
    rev_nonoverlap = (rev["fl3"] - fwd["fl5"]) if (fwd and rev) else 0
    nonoverlap_ok  = (fwd_nonoverlap >= 8 and rev_nonoverlap >= 8)

    return dict(
        tmpl_len=len(T), pos=pos, orig=O, mut=M,
        n_mm=n_mm, fwd=fwd, rev=rev,
        nonoverlap_ok=nonoverlap_ok,
        fwd_nonoverlap=fwd_nonoverlap,
        rev_nonoverlap=rev_nonoverlap,
        codon_info=None,
    )


# ----------------------- combined design ----------------------------

def _annotate_fwd(seq, muts_sorted, primer_start):
    result = []
    seq_i = 0
    tmpl_pos = primer_start
    for ms, me, orig, new in muts_sorted:
        pre = ms - tmpl_pos
        result.append(seq[seq_i:seq_i + pre].lower())
        seq_i += pre
        if new:
            result.append("[" + seq[seq_i:seq_i + len(new)] + "]")
            seq_i += len(new)
        tmpl_pos = me
    result.append(seq[seq_i:].lower())
    return "".join(result)


def _annotate_rev(seq, muts_sorted, fl5_rev):
    result = []
    seq_i = 0
    n = len(muts_sorted)

    result.append(seq[seq_i:seq_i + fl5_rev].lower())
    seq_i += fl5_rev

    for k in range(n - 1, -1, -1):
        ms, me, orig, new = muts_sorted[k]
        result.append("[" + seq[seq_i:seq_i + len(new)] + "]")
        seq_i += len(new)
        if k > 0:
            inter = muts_sorted[k][0] - muts_sorted[k - 1][1]
            if inter > 0:
                result.append(seq[seq_i:seq_i + inter].lower())
                seq_i += inter

    result.append(seq[seq_i:].lower())
    return "".join(result)


def design_combined(template, mutations):
    T = validate_dna(template)

    muts = []
    for pos, orig, new in mutations:
        O = validate_dna(orig) if orig else ""
        M = validate_dna(new)  if new  else ""
        ms, me = pos - 1, pos - 1 + len(O)
        if me > len(T):
            raise ValueError("Position %d + length %d exceeds template" % (pos, len(O)))
        if O and T[ms:me] != O:
            raise ValueError("Position %d: template '%s' != expected '%s'"
                             % (pos, T[ms:me], O))
        muts.append((ms, me, O, M))

    muts.sort(key=lambda x: x[0])

    for i in range(len(muts) - 1):
        if muts[i][1] > muts[i + 1][0]:
            raise ValueError(
                "Mutations at nt %d and nt %d overlap -- cannot combine"
                % (muts[i][0] + 1, muts[i + 1][0] + 1)
            )

    first_ms = muts[0][0]
    last_me  = muts[-1][1]
    n_mm = sum(count_mm(O, M) for _, _, O, M in muts)

    def build_seq(fl5, fl3, fwd):
        s, e = (first_ms - fl5, last_me + fl3) if fwd else (first_ms - fl3, last_me + fl5)
        if s < 0 or e > len(T):
            return None, None, None
        region = T[s:muts[0][0]]
        for i, (ms, me, orig, new) in enumerate(muts):
            region += new
            if i < len(muts) - 1:
                region += T[me:muts[i + 1][0]]
        region += T[muts[-1][1]:e]
        return (region if fwd else revcomp(region)), s, e

    def get_candidates(fwd):
        ML5, ML3 = 4, 8
        lim5 = min(first_ms if fwd else (len(T) - last_me), 70)
        lim3 = min((len(T) - last_me) if fwd else first_ms, 70)
        cands = []
        for fl5 in range(ML5, lim5 + 1):
            for fl3 in range(ML3, lim3 + 1):
                seq, s, e = build_seq(fl5, fl3, fwd)
                if seq is None or seq[-1] not in "GC":
                    continue
                tm = calc_tm(seq, n_mm)
                if tm >= 78.0:
                    cands.append({"seq": seq, "tm": round(tm, 2), "len": len(seq),
                                  "gc": round(gc_pct(seq), 1),
                                  "fl5": fl5, "fl3": fl3, "ts": s + 1, "te": e})
        cands.sort(key=lambda x: x["len"])
        return cands

    fwd_list = get_candidates(fwd=True)
    rev_list = get_candidates(fwd=False)

    best_fwd = best_rev = None
    best_total = float("inf")
    for f in fwd_list:
        for r in rev_list:
            if f["fl3"] < r["fl5"] + 8: continue
            if r["fl3"] < f["fl5"] + 8: continue
            t = f["len"] + r["len"]
            if t < best_total:
                best_total, best_fwd, best_rev = t, f, r

    fn = (best_fwd["fl3"] - best_rev["fl5"]) if (best_fwd and best_rev) else 0
    rn = (best_rev["fl3"] - best_fwd["fl5"]) if (best_fwd and best_rev) else 0

    return dict(
        tmpl_len=len(T), muts=muts,
        first_pos=first_ms + 1, last_pos=last_me, span=last_me - first_ms,
        n_mm=n_mm, fwd=best_fwd, rev=best_rev,
        nonoverlap_ok=(fn >= 8 and rn >= 8),
        fwd_nonoverlap=fn, rev_nonoverlap=rn,
    )


# --------------------------- primer naming ---------------------------

def _mut_label_single(pos, orig, new, codon_info=None):
    """Return the mutation part of a primer name (no gene prefix, no F/R suffix)."""
    if codon_info:
        return "%s%d%s" % (codon_info["orig_aa"], codon_info["aa_pos"], codon_info["new_aa"])
    if not orig:
        return "nt%dins%s" % (pos, new)
    if not new:
        return "nt%ddel%s" % (pos, orig)
    if len(orig) == 1 and len(new) == 1:
        return "%s%d%s" % (orig, pos, new)
    return "nt%d%sto%s" % (pos, orig, new)


def _mut_label_combined(muts_sorted, codon_infos=None):
    """Return underscore-joined mutation labels for a combined primer pair."""
    parts = []
    for i, (ms, me, orig, new) in enumerate(muts_sorted):
        ci = codon_infos[i] if codon_infos else None
        parts.append(_mut_label_single(ms + 1, orig, new, ci))
    return "_".join(parts)


def primer_name(gene, label, is_fwd):
    """Assemble final primer name: [Gene_]{label}_F/R"""
    direction = "F" if is_fwd else "R"
    if gene:
        return "%s_%s_%s" % (gene, label, direction)
    return "%s_%s" % (label, direction)


# --------------------------- output ----------------------------------

def show(r, gene_name=""):
    if r["orig"] and r["mut"]:
        kind = "SUBSTITUTION"
    elif not r["orig"]:
        kind = "INSERTION"
    else:
        kind = "DELETION"

    W = 68
    sep  = "=" * W
    thin = "-" * W

    ci = r.get("codon_info")
    label = _mut_label_single(r["pos"], r["orig"], r["mut"], ci)
    fwd_name = primer_name(gene_name, label, is_fwd=True)
    rev_name = primer_name(gene_name, label, is_fwd=False)

    print()
    print(sep)
    print("  SDM DOUBLE PRIMER DESIGN RESULTS")
    print(sep)
    print("  Template  : %d bp" % r["tmpl_len"])
    print("  Mutation  : %s  at nt position %d" % (kind, r["pos"]))
    print("  Original  : %s" % (r["orig"] or "(none)"))
    print("  New seq   : %s" % (r["mut"]  or "(none)"))
    print("  Mismatches: %d bp  (%%mismatch = %d/N x 100)" % (r["n_mm"], r["n_mm"]))

    if ci:
        print()
        print(thin)
        print("  CODON SELECTION  (E. coli K-12 MG1655)")
        print(thin)
        print("  AA position : %d" % ci["aa_pos"])
        print("  Substitution: %s (%s)  -->  %s (%s)"
              % (ci["orig_aa"], AA_NAMES.get(ci["orig_aa"], ""),
                 ci["new_aa"],  AA_NAMES.get(ci["new_aa"],  "")))
        print("  Codon change: %s  -->  %s" % (ci["orig_codon"], ci["new_codon"]))
        print("  CDS start   : nt %d  (codon at nt %d-%d)"
              % (ci["cds_start"], ci["nt_pos"], ci["nt_pos"] + 2))
        print()
        print("  All codons for %s (%s) by E. coli K-12 frequency:"
              % (ci["new_aa"], AA_NAMES.get(ci["new_aa"], "")))
        for codon, freq in ci["ranked"]:
            tag = "  <-- selected" if codon == ci["new_codon"] else ""
            print("    %s   %5.1f / 1000%s" % (codon, freq, tag))

    for key, dir_label, is_fwd, pname in [
        ("fwd", "FORWARD PRIMER   5'->3'  (sense strand)",    True,  fwd_name),
        ("rev", "REVERSE PRIMER   5'->3'  (antisense strand)", False, rev_name),
    ]:
        p = r[key]
        print()
        print(thin)
        print("  %s" % dir_label)
        print("  Primer name : %s" % pname)
        print(thin)
        if p is None:
            print("  [NOT FOUND]  No primer satisfying all criteria within max length.")
            print("  Tip: provide more flanking template sequence (>= 30 nt each side).")
            continue
        print("  Sequence : 5'-%s-3'" % p["seq"])
        print("  Length   : %d bp     Tm: %.2f degC     %%GC: %.1f%%"
              % (p["len"], p["tm"], p["gc"]))
        print("  5' flank : %d nt  (mutation %d nt from 5' end)" % (p["fl5"], p["fl5"]))
        print("  3' flank : %d nt  (mutation %d nt from 3' end)" % (p["fl3"], p["fl3"]))
        print("  Template : nt %d - %d" % (p["ts"], p["te"]))

    print()
    print(thin)
    print("  VALIDATION CHECKS")
    print(thin)

    f, v = r["fwd"], r["rev"]
    fn = r["fwd_nonoverlap"]
    rn = r["rev_nonoverlap"]
    checks = [
        ("Mutation present in both primers",
             f is not None and v is not None),
        ("Fwd: mutation >= 4 nt from 5' end",
             f is not None and f["fl5"] >= 4),
        ("Rev: mutation >= 4 nt from 5' end",
             v is not None and v["fl5"] >= 4),
        ("Fwd: mutation >= 8 nt from 3' end",
             f is not None and f["fl3"] >= 8),
        ("Rev: mutation >= 8 nt from 3' end",
             v is not None and v["fl3"] >= 8),
        ("Fwd: 3' terminal base is G or C",
             f is not None and f["seq"][-1] in "GC"),
        ("Rev: 3' terminal base is G or C",
             v is not None and v["seq"][-1] in "GC"),
        ("Fwd: Tm >= 78 degC",
             f is not None and f["tm"] >= 78.0),
        ("Rev: Tm >= 78 degC",
             v is not None and v["tm"] >= 78.0),
        ("Fwd 3' non-overlapping = fl3_fwd(%d) - fl5_rev(%d) = %d >= 8"
             % (f["fl3"] if f else 0, v["fl5"] if v else 0, fn),
             fn >= 8),
        ("Rev 3' non-overlapping = fl3_rev(%d) - fl5_fwd(%d) = %d >= 8"
             % (v["fl3"] if v else 0, f["fl5"] if f else 0, rn),
             rn >= 8),
    ]
    all_pass = True
    for desc, ok in checks:
        print("  %s  %s" % ("[OK]  " if ok else "[FAIL]", desc))
        if not ok:
            all_pass = False

    print()
    print("  >> All criteria satisfied." if all_pass
          else "  >> One or more criteria NOT met. Check [FAIL] items above.")
    print(sep)
    print()


def show_combined(r, label="COMBINED PRIMER DESIGN", codon_infos=None, gene_name=""):
    W = 68
    sep  = "=" * W
    thin = "-" * W
    muts = r["muts"]

    combined_label = _mut_label_combined(muts, codon_infos)
    fwd_name = primer_name(gene_name, combined_label, is_fwd=True)
    rev_name = primer_name(gene_name, combined_label, is_fwd=False)

    print()
    print(sep)
    print("  %s  (%d mutations)" % (label, len(muts)))
    print(sep)
    print("  Template  : %d bp" % r["tmpl_len"])
    print("  Mutations span : nt %d - %d  (%d bp)"
          % (r["first_pos"], r["last_pos"], r["span"]))
    print("  Total mismatches: %d" % r["n_mm"])
    print()
    for i, (ms, me, orig, new) in enumerate(muts):
        ci = codon_infos[i] if codon_infos else None
        aa_note = ""
        if ci:
            aa_note = "  [AA%d %s->%s, codon %s->%s (E.coli K-12)]" % (
                ci["aa_pos"], ci["orig_aa"], ci["new_aa"],
                ci["orig_codon"], ci["new_codon"])
        kind = ("SUB" if orig and new else "INS" if not orig else "DEL")
        print("    [%d] nt %d  %s  %s -> %s%s"
              % (i + 1, ms + 1, kind, orig or "(none)", new or "(none)", aa_note))

    for key, dir_label, fwd_flag, pname in [
        ("fwd", "FORWARD 5'->3' (sense)",    True,  fwd_name),
        ("rev", "REVERSE 5'->3' (antisense)", False, rev_name),
    ]:
        p = r[key]
        print()
        print(thin)
        print("  %s" % dir_label)
        print("  Primer name : %s" % pname)
        print(thin)
        if p is None:
            print("  [NOT FOUND]  Mutations may be too far apart or near template ends.")
            print("  Tip: extend flanking template sequence or design separate primers.")
            continue
        if fwd_flag:
            ann = _annotate_fwd(p["seq"], muts, p["ts"] - 1)
        else:
            ann = _annotate_rev(p["seq"], muts, p["fl5"])
        print("  Annotated: 5'-%s-3'" % ann)
        print("  Sequence : 5'-%s-3'" % p["seq"])
        print("  Length   : %d bp     Tm: %.2f degC     %%GC: %.1f%%"
              % (p["len"], p["tm"], p["gc"]))
        print("  5' flank : %d nt  |  3' flank : %d nt" % (p["fl5"], p["fl3"]))
        print("  Template : nt %d - %d" % (p["ts"], p["te"]))

    print()
    print(thin)
    print("  VALIDATION")
    print(thin)
    f, v = r["fwd"], r["rev"]
    fn, rn = r["fwd_nonoverlap"], r["rev_nonoverlap"]
    checks = [
        ("All mutations in both primers",       f is not None and v is not None),
        ("Fwd: first mut >= 4 nt from 5' end",  f is not None and f["fl5"] >= 4),
        ("Rev: first mut >= 4 nt from 5' end",  v is not None and v["fl5"] >= 4),
        ("Fwd: last mut >= 8 nt from 3' end",   f is not None and f["fl3"] >= 8),
        ("Rev: last mut >= 8 nt from 3' end",   v is not None and v["fl3"] >= 8),
        ("Fwd: 3' terminal is G or C",          f is not None and f["seq"][-1] in "GC"),
        ("Rev: 3' terminal is G or C",          v is not None and v["seq"][-1] in "GC"),
        ("Fwd: Tm >= 78 degC",                  f is not None and f["tm"] >= 78.0),
        ("Rev: Tm >= 78 degC",                  v is not None and v["tm"] >= 78.0),
        ("Fwd 3' non-overlapping = fl3(%d)-fl5_rev(%d) = %d >= 8"
             % (f["fl3"] if f else 0, v["fl5"] if v else 0, fn), fn >= 8),
        ("Rev 3' non-overlapping = fl3(%d)-fl5_fwd(%d) = %d >= 8"
             % (v["fl3"] if v else 0, f["fl5"] if f else 0, rn), rn >= 8),
    ]
    all_pass = True
    for desc, ok in checks:
        print("  %s  %s" % ("[OK]  " if ok else "[FAIL]", desc))
        if not ok:
            all_pass = False
    print()
    print("  >> All criteria satisfied." if all_pass
          else "  >> One or more criteria NOT met.")
    print(sep)
    print()


# --------------------------- input helpers ---------------------------

def parse_dna_mut(s):
    """'<pos> <orig|-> <new|->'  ->  (pos, orig, new)"""
    parts = s.strip().split()
    if len(parts) != 3:
        raise ValueError("Need exactly 3 fields: position  original  new")
    pos = int(parts[0])
    if pos < 1:
        raise ValueError("Position must be >= 1")
    orig = "" if parts[1] == "-" else parts[1].upper()
    new  = "" if parts[2] == "-" else parts[2].upper()
    if not orig and not new:
        raise ValueError("original and new cannot both be '-'")
    return pos, orig, new


def parse_aa_mut(s):
    """'<aa_pos> <orig_aa> <new_aa>'  ->  (aa_pos, orig_aa, new_aa)"""
    parts = s.strip().split()
    if len(parts) != 3:
        raise ValueError("Need exactly 3 fields: aa_position  original_AA  new_AA")
    pos = int(parts[0])
    if pos < 1:
        raise ValueError("AA position must be >= 1")
    orig = parts[1].upper()
    new  = parts[2].upper()
    return pos, orig, new


# --------------------------- SSM -------------------------------------

_ALL_AAS = "ACDEFGHIKLMNPQRSTVWY"


def show_ssm(results, gene_name, wt_aa, aa_pos, wt_codon, tmpl_len):
    W    = 68
    sep  = "=" * W
    thin = "-" * W

    pass_count  = sum(1 for e in results if e["success"] and e["pass"])
    check_count = sum(1 for e in results if e["success"] and not e["pass"])
    err_count   = sum(1 for e in results if not e["success"])

    print()
    print(sep)
    print("  SITE SATURATION MUTAGENESIS")
    print(sep)
    print("  Site     : AA%d  %s (%s)  wild-type codon: %s"
          % (aa_pos, wt_aa, AA_NAMES.get(wt_aa, "?"), wt_codon))
    print("  Template : %d bp" % tmpl_len)
    print("  Results  : %d PASS  /  %d CHECK  /  %d ERROR  (of 19)"
          % (pass_count, check_count, err_count))
    print(sep)
    print()

    # Dynamic column width for primer names
    name_w = max(
        (len(e["fwd_name"]) for e in results if e.get("success")),
        default=18,
    )
    name_w = max(name_w, 18)

    hdr_fmt = "  %-10s %-6s %-3s  %-*s  %5s  %7s  %s"
    row_fmt  = "  %-10s %-6s %-3s  %-*s  %5s  %7s  %s"
    row_fmt2 = "  %-10s %-6s %-3s  %-*s  %5s  %7s"

    print(hdr_fmt % ("AA", "Codon", "Dir", name_w, "Primer Name", "Len", "Tm(°C)", "OK?"))
    print("  " + thin[2:])

    for entry in results:
        aa_label = "%s (%s)" % (entry["target_aa"], AA_NAMES.get(entry["target_aa"], "?"))

        if not entry["success"]:
            print("  %-10s %-6s  ERROR: %s" % (aa_label, "-", entry["error"]))
            print()
            continue

        ci     = entry["codon_info"]
        r      = entry["r"]
        status = "PASS" if entry["pass"] else "CHECK"

        fwd = r["fwd"]
        if fwd:
            print(row_fmt % (aa_label, ci["new_codon"], "F",
                             name_w, entry["fwd_name"],
                             fwd["len"], "%.1f" % fwd["tm"], status))
        else:
            print(row_fmt % (aa_label, ci["new_codon"], "F",
                             name_w, "[NOT FOUND]", "-", "-", "CHECK"))

        rev = r["rev"]
        if rev:
            print(row_fmt2 % ("", "", "R",
                              name_w, entry["rev_name"],
                              rev["len"], "%.1f" % rev["tm"]))
        else:
            print(row_fmt2 % ("", "", "R", name_w, "[NOT FOUND]", "-", "-"))

    print()
    print("  " + thin[2:])
    print("  PASS: %d  |  CHECK: %d  |  ERROR: %d"
          % (pass_count, check_count, err_count))
    print(sep)
    if check_count or err_count:
        print()
        print("  Tip: CHECK rows need attention. Provide >= 30 nt flanking sequence")
        print("       on each side of the mutation site and re-run.")
    print()


def run_ssm(tmpl, cds_start, gene_name=""):
    print()
    print("  Site Saturation Mutagenesis")
    print("  Enter the amino acid position to saturate (1-based).")
    print("  The wild-type residue is auto-detected from the template.")
    print()

    while True:
        try:
            aa_pos = int(input("  AA position > ").strip())
            if aa_pos < 1:
                raise ValueError("Must be >= 1")
            nt_off = cds_start + (aa_pos - 1) * 3 - 1
            if nt_off < 0 or nt_off + 3 > len(tmpl):
                raise ValueError(
                    "AA position %d (nt %d-%d) is outside the template (%d bp)"
                    % (aa_pos, nt_off + 1, nt_off + 3, len(tmpl))
                )
            wt_codon = tmpl[nt_off:nt_off + 3]
            wt_aa    = translate_codon(wt_codon)
            if wt_aa == "*":
                raise ValueError("Position %d encodes a stop codon (%s)" % (aa_pos, wt_codon))
            if wt_aa == "?":
                raise ValueError("Codon '%s' at position %d is not standard" % (wt_codon, aa_pos))
            print("  Wild-type at AA%d: %s = %s (%s)"
                  % (aa_pos, wt_codon, wt_aa, AA_NAMES.get(wt_aa, "?")))
            break
        except ValueError as e:
            print("  Error: %s -- try again." % e)

    nt_pos = cds_start + (aa_pos - 1) * 3   # 1-based nt position of codon start

    print()
    print("  Designing primers for all 19 substitutions at %s%d ..." % (wt_aa, aa_pos))

    results = []
    for target_aa in _ALL_AAS:
        if target_aa == wt_aa:
            continue
        try:
            new_codon   = best_ecoli_codon(target_aa)
            codon_info  = {
                "aa_pos":     aa_pos,
                "orig_aa":    wt_aa,
                "new_aa":     target_aa,
                "orig_codon": wt_codon,
                "new_codon":  new_codon,
                "ranked":     ecoli_codons_ranked(target_aa),
                "nt_pos":     nt_pos,
                "cds_start":  cds_start,
            }
            r = design(tmpl, nt_pos, wt_codon, new_codon)
            r["codon_info"] = codon_info

            mut_lbl  = "%s%d%s" % (wt_aa, aa_pos, target_aa)
            fwd_name = primer_name(gene_name, mut_lbl, is_fwd=True)
            rev_name = primer_name(gene_name, mut_lbl, is_fwd=False)
            ok = bool(r["fwd"] and r["rev"] and r["nonoverlap_ok"]
                      and r["fwd"]["tm"] >= 78.0 and r["rev"]["tm"] >= 78.0)
            results.append(dict(
                target_aa=target_aa, codon_info=codon_info,
                r=r, fwd_name=fwd_name, rev_name=rev_name,
                success=True, **{"pass": ok},
            ))
        except Exception as e:
            results.append(dict(target_aa=target_aa, error=str(e), success=False))

    show_ssm(results, gene_name, wt_aa, aa_pos, wt_codon, len(tmpl))

    # Optional detailed view
    while True:
        print("Show full primer details? [all / <AA letter> / n]: ", end="")
        ans = input().strip().lower()
        if ans in ("", "n", "no"):
            break
        if ans == "all":
            for entry in results:
                if entry["success"]:
                    print()
                    print("  === %s%d%s ===" % (wt_aa, aa_pos, entry["target_aa"]))
                    show(entry["r"], gene_name=gene_name)
            break
        if ans.upper() in _ALL_AAS:
            target = ans.upper()
            found = next((e for e in results if e["target_aa"] == target), None)
            if found is None:
                print("  '%s' is the wild-type or not a standard amino acid." % target.upper())
            elif not found["success"]:
                print("  Error for %s%d%s: %s" % (wt_aa, aa_pos, target, found["error"]))
            else:
                show(found["r"], gene_name=gene_name)
        else:
            print("  Enter 'all', a 1-letter AA code (e.g. 'A'), or 'n' to skip.")


# --------------------------- main ------------------------------------

def _resolve_mutation(tmpl, raw, cds_start):
    parts = raw.strip().split()
    if not parts:
        raise ValueError("Empty input")

    if parts[0].lower() == "aa":
        if len(parts) != 4:
            raise ValueError("AA format: aa <position> <orig_AA> <new_AA>")
        aa_pos = int(parts[1])
        orig_aa, new_aa = parts[2].upper(), parts[3].upper()
        pos, orig, new, ci = aa_mutation_to_dna(tmpl, aa_pos, orig_aa, new_aa, cds_start)
        return pos, orig, new, ci
    else:
        if parts[0].lower() == "nt":
            parts = parts[1:]
        pos, orig, new = parse_dna_mut(" ".join(parts))
        return pos, orig, new, None


def run_batch(tmpl, cds_start, gene_name=""):
    print()
    print("  Enter mutations one per line, blank line when done.")
    print("  Formats:")
    print("    <pos> <orig> <new>          DNA level  e.g.  73 TAC GCG")
    print("    nt <pos> <orig> <new>       DNA level  e.g.  nt 73 TAC GCG")
    print("    aa <pos> <orig_AA> <new_AA> AA level   e.g.  aa 25 Y A")
    print("    (use - for empty orig/new in insertions/deletions)")
    print()

    collected = []
    idx = 1
    while True:
        raw = input("  Mutation %d > " % idx).strip()
        if raw == "":
            if not collected:
                print("  No mutations entered.")
                continue
            break
        try:
            pos, orig, new, ci = _resolve_mutation(tmpl, raw, cds_start)
            collected.append((pos, orig, new, ci))
            print("    -> nt %d  %s -> %s  %s"
                  % (pos, orig or "(ins)", new or "(del)",
                     ("[AA%d %s->%s codon %s->%s]" % (ci["aa_pos"], ci["orig_aa"],
                      ci["new_aa"], ci["orig_codon"], ci["new_codon"])) if ci else ""))
            idx += 1
        except (ValueError, IndexError) as e:
            print("  Error: %s  -- try again." % e)

    n = len(collected)
    print()
    print("  %d mutation%s collected." % (n, "s" if n > 1 else ""))

    print()
    print("=" * 68)
    print("  INDIVIDUAL PRIMER PAIRS  (%d)" % n)
    print("=" * 68)
    for i, (pos, orig, new, ci) in enumerate(collected):
        print()
        print("  --- Mutation %d of %d: nt %d  %s -> %s ---"
              % (i + 1, n, pos, orig or "(ins)", new or "(del)"))
        try:
            r = design(tmpl, pos, orig, new)
            r["codon_info"] = ci
            show(r, gene_name=gene_name)
        except ValueError as e:
            print("  [ERROR] %s\n" % e)

    if n < 2:
        return

    print()
    print("=" * 68)
    print("  COMBINED PRIMER PAIR  (all %d mutations in one primer set)" % n)
    print("=" * 68)
    try:
        mut_list = [(pos, orig, new) for pos, orig, new, _ in collected]
        codon_infos = [ci for _, _, _, ci in collected]
        rc = design_combined(tmpl, mut_list)
        show_combined(rc, label="COMBINED PRIMER DESIGN",
                      codon_infos=codon_infos, gene_name=gene_name)
    except ValueError as e:
        print()
        print("  [COMBINED NOT POSSIBLE] %s" % e)
        print()


def run_once(tmpl, cds_start, gene_name=""):
    print()
    print("Mutation mode:")
    print("  [1] DNA level   -- specify nucleotide position and sequence")
    print("  [2] AA level    -- specify amino acid position (E. coli codon auto-selected)")

    while True:
        mode = input("\n  Choose [1/2]: ").strip()
        if mode in ("1", "2"):
            break
        print("  Please enter 1 or 2.")

    codon_info = None

    if mode == "1":
        print()
        print("Format:  <1-based nt position>  <original>  <new>  (use - for empty)")
        print("  Examples:")
        print("    152 A G         A -> G at position 152")
        print("    80  ATG CTG     ATG -> CTG starting at position 80")
        print("    45  -  GAATTC   insert GAATTC before position 45")
        print("    30  TGCA -      delete TGCA starting at position 30")

        while True:
            try:
                pos, orig, new = parse_dna_mut(input("\n  > "))
                break
            except (ValueError, IndexError) as e:
                print("  Error: %s  -- try again." % e)

    else:
        print()
        print("Format:  <AA position>  <original AA>  <new AA>  (1-letter codes)")
        print("  Example:  52 A V    (Ala52Val)")

        while True:
            try:
                aa_pos, orig_aa, new_aa = parse_aa_mut(input("\n  > "))
                pos, orig, new, codon_info = aa_mutation_to_dna(
                    tmpl, aa_pos, orig_aa, new_aa, cds_start
                )
                break
            except (ValueError, IndexError) as e:
                print("  Error: %s  -- try again." % e)

    try:
        result = design(tmpl, pos, orig, new)
        result["codon_info"] = codon_info
        show(result, gene_name=gene_name)
    except ValueError as e:
        print("\n  Error: %s\n" % e)


def main():
    print()
    print("  +--------------------------------------------------+")
    print("  |   SDM Double Primer Designer  v3                 |")
    print("  |   Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch    |")
    print("  |   AA mode: E. coli K-12 codon optimization      |")
    print("  +--------------------------------------------------+")

    # -- template --
    while True:
        print()
        print("Template DNA sequence (5'->3', paste and press Enter):")
        tmpl = input("  > ").strip()
        try:
            clean = validate_dna(tmpl)
            if len(clean) < 20:
                print("  Template too short; please provide >= 20 bp.")
                continue
            break
        except ValueError as e:
            print("  Error: %s" % e)

    tmpl = clean
    print("  Template accepted: %d bp" % len(tmpl))

    # -- gene name --
    print()
    print("Gene name (used for primer naming, e.g. MaEgtB; leave blank to skip):")
    gene_name = input("  Gene name []: ").strip()
    if gene_name:
        print("  Gene: %s  (primers will be named %s_<mut>_F / %s_<mut>_R)"
              % (gene_name, gene_name, gene_name))
    else:
        print("  No gene name set.")

    # -- CDS start --
    print()
    print("CDS start position (1-based nt, for AA mode; default = 1):")
    while True:
        raw = input("  CDS start [1]: ").strip()
        if raw == "":
            cds_start = 1
            break
        try:
            cds_start = int(raw)
            if cds_start >= 1:
                break
            print("  Must be >= 1.")
        except ValueError:
            print("  Enter a number.")

    # -- design mode --
    print()
    print("Design mode:")
    print("  [1] Single mutation")
    print("  [2] Multiple mutations (individual + combined primer sets)")
    print("  [3] Site Saturation Mutagenesis -- all 19 AA substitutions at one position")

    while True:
        mode = input("\n  Choose [1/2/3]: ").strip()
        if mode in ("1", "2", "3"):
            break
        print("  Please enter 1, 2, or 3.")

    if mode == "1":
        run_once(tmpl, cds_start, gene_name=gene_name)
        while True:
            ans = input("Design another mutation on the same template? [y/N]: ").strip().lower()
            if ans not in ("y", "yes"):
                break
            run_once(tmpl, cds_start, gene_name=gene_name)
    elif mode == "2":
        run_batch(tmpl, cds_start, gene_name=gene_name)
        while True:
            ans = input("Design another batch on the same template? [y/N]: ").strip().lower()
            if ans not in ("y", "yes"):
                break
            run_batch(tmpl, cds_start, gene_name=gene_name)
    else:
        run_ssm(tmpl, cds_start, gene_name=gene_name)
        while True:
            ans = input("Design SSM at another position on the same template? [y/N]: ").strip().lower()
            if ans not in ("y", "yes"):
                break
            run_ssm(tmpl, cds_start, gene_name=gene_name)


if __name__ == "__main__":
    main()
