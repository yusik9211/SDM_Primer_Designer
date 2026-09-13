# SDM Primer Designer

An automated double primer design tool for Site-Directed Mutagenesis (SDM).  
Given a template DNA sequence and desired mutations, it automatically designs forward/reverse primer pairs that satisfy all NEB-style SDM design rules.

---

## Files

| File | Description |
|------|-------------|
| `sdm_primer_designer.html` | **Web version (recommended)** — open in browser, no installation; Excel/TXT download |
| `sdm_primer_designer.py`   | **Python version** — run in terminal, requires Python 3 |

Previous versions (v1, v2) are available in the [`legacy/`](legacy/) folder.

---

## Quick Start

### HTML version
Open `sdm_primer_designer_v3.html` in your browser — no installation required.  
(Internet connection required for Excel export via SheetJS CDN)

### Python version
```bash
python sdm_primer_designer_v3.py
```
Requires Python 3.6+. No external packages needed.

---

## Usage — HTML Version

### Step 1 — Enter Template DNA Sequence
- Paste your template DNA sequence in the 5'→3' direction (A/T/G/C/N only)
- Enter a **Gene name** (e.g. `MaEgtB`) — automatically reflected in primer names
- Set the CDS start position (required for AA mutation mode; default = 1)

### Step 2 — Add Mutations
Select a mutation type and add entries. Multiple mutations can be entered at once.

#### DNA Mutation
| Field | Description | Example |
|-------|-------------|---------|
| Position | Mutation start position (1-based) | `298` |
| Original | Original sequence (`-` = insertion) | `GGC` |
| New | New sequence (`-` = deletion) | `GCC` |

#### Amino Acid (AA) Mutation
| Field | Description | Example |
|-------|-------------|---------|
| AA Position | Amino acid position (1-based) | `100` |
| Orig AA | Original amino acid (1-letter code) | `G` |
| New AA | New amino acid (1-letter code) | `A` |

> In AA mutation mode, the optimal codon is automatically selected based on *E. coli* K-12 MG1655 codon usage frequency.

#### Upload a Mutant List (Excel/CSV)
Instead of adding rows one by one, you can upload a spreadsheet listing all your desired mutants at once:
- Click **Upload Mutant List (Excel/CSV)** in Step 2 and select an `.xlsx`, `.xls`, or `.csv` file
- Each mutant is written as **OrigAA + Position + NewAA** notation in a single cell, e.g. `A30G` (Ala30→Gly)
- Any column/row layout works — every matching cell in the sheet is detected automatically
- Each recognized mutation is added as an independent AA mutation row (its own primer pair; rows are **not** combined into one multi-mutation construct)
- Click **Download template** for a ready-to-fill example file

### Step 3 — Click "Design Primers"
- **Individual Primer Pairs**: designs an independent primer pair for each mutation
- **Combined Primer Pair**: incorporates all mutations into a single primer pair (requires ≥ 2 mutations)

### Results and Download
- Each primer is labeled with a name tag (e.g. `MaEgtB_G100A_F`)
- A Validation Checklist confirms all design rules are satisfied
- **Download TXT**: saves all results as a plain text file
- **Download Excel**: exports primer name, sequence, Tm, and %GC to an Excel file

---

## Usage — Python Version

```
python sdm_primer_designer_v3.py
```

Follow the prompts in order:

```
Template DNA sequence: ATGCGT...    ← paste DNA sequence
Gene name []:          MaEgtB       ← gene name (press Enter to skip)
CDS start [1]:         1            ← CDS start position

Design mode:
  [1] Single mutation
  [2] Multiple mutations
Choose [1/2]: 2
```

#### Mutation Input Format
```
# DNA mutation
<position> <original> <new>
73 TAC GCG          → nt73 TAC→GCG (substitution)
45 - GAATTC         → insert GAATTC at nt45
30 TGCA -           → delete TGCA starting at nt30

# DNA mutation (nt prefix also accepted)
nt 73 TAC GCG

# Amino acid mutation
aa <aa_position> <orig_AA> <new_AA>
aa 100 G A          → Gly100 → Ala
aa 25 Y A           → Tyr25 → Ala
```

---

## Primer Naming Convention

```
{Gene}_{Mutation}_{Direction}
```

| Mutation type | Input example | Forward name | Reverse name |
|---------------|--------------|--------------|--------------|
| AA substitution | `aa 100 G A` | `MaEgtB_G100A_F` | `MaEgtB_G100A_R` |
| Single nt substitution | `300 T A` | `MaEgtB_T300A_F` | `MaEgtB_T300A_R` |
| Multi-nt substitution | `80 ATG CTG` | `MaEgtB_nt80ATGtoCTG_F` | `MaEgtB_nt80ATGtoCTG_R` |
| Insertion | `45 - GAATTC` | `MaEgtB_nt45insGAATTC_F` | `MaEgtB_nt45insGAATTC_R` |
| Deletion | `30 TGCA -` | `MaEgtB_nt30delTGCA_F` | `MaEgtB_nt30delTGCA_R` |
| Multiple (combined) | G100A + L200P | `MaEgtB_G100A_L200P_F` | `MaEgtB_G100A_L200P_R` |

If no gene name is provided, the format is shortened (e.g. `G100A_F`).

---

## Primer Design Rules

All designed primers simultaneously satisfy the following seven criteria:

| Rule | Constraint |
|------|-----------|
| Mutation in both primers | Both forward and reverse primers span the mutation site |
| 5′ flank | Mutation site ≥ 4 nt from the 5′ terminus |
| 3′ flank | Mutation site ≥ 8 nt from the 3′ terminus |
| Non-overlapping 3′ extension | ≥ 8 nt at the 3′ end that do not overlap the partner primer's 3′ region |
| 3′ terminal base | Must be G or C |
| Melting temperature | Tm ≥ 78 °C |
| Tm formula | `Tm = 81.5 + 0.41(%GC) − 675/N − %mismatch` |

The shortest primer pair satisfying all constraints is automatically selected.

---

## *E. coli* Codon Optimization

In AA mutation mode, the codon with the highest usage frequency (per 1,000 codons) among all synonymous codons for the target amino acid is automatically chosen.

- Source: Codon Usage Database (Kazusa), *E. coli* K-12 MG1655
- In the HTML version, click the **Codon selection** toggle in the results panel to view all candidate codons and their frequencies.

---

## Output Example

```
  FORWARD PRIMER   5'->3'  (sense strand)
  Primer name : MaEgtB_G100A_F
  ----------------------------------------------------------------
  Sequence : 5'-GCCTGATCGTtGCCATGCGT-3'
  Length   : 20 bp     Tm: 79.34 degC     %GC: 55.0%
  5' flank : 10 nt  |  3' flank : 8 nt
  Template : nt 289 - 308
```

---

## Notes

- The template sequence should have **at least 30 nt of flanking sequence on each side** of the mutation site. Insufficient flanking may prevent a valid primer from being found.
- In Combined mode with multiple mutations, overlapping mutation sites will cause design failure.
- Excel download in the HTML version requires an internet connection (SheetJS CDN).

---

## Reference

Zheng, L., U. Baumann, and Jean-Louis Reymond. (2004).
An efficient one-step site-directed and site-saturation mutagenesis protocol.
*Nucleic Acids Res.* 32(14): e115.

---

## Copyright

Copyright &copy; 2026 Yusik Kim, Korea University  
Contact: [yusik9211@korea.ac.kr](mailto:yusik9211@korea.ac.kr)
