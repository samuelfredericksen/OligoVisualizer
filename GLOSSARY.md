# Oligo viewer

Vocabulary for describing synthetic oligonucleotides as IDT writes and sells them.

## Chain

**Oligonucleotide**:
A short synthetic nucleic acid chain, written 5' end first.
_Avoid_: Oligo (fine in speech, not in labels), primer, probe

**Nucleotide**:
One position in the chain: a nucleobase plus its sugar and its linkage to the next position.
_Avoid_: Base (a base is only part of a nucleotide), residue

**Nucleobase**:
The A, C, G, T, U or inosine ring at a position; the thing that gets the colour.
_Avoid_: Base letter, character

**Linkage**:
The bond joining one nucleotide to the next. Two kinds occur here: phosphate and
phosphorothioate.
_Avoid_: Bond, backbone bond

**Phosphorothioate linkage**:
A phosphate linkage with one non-bridging oxygen replaced by sulfur; written `*` after the
nucleotide, it blocks nuclease digestion.
_Avoid_: Thiolated linkage, PS bond, thioate

**Phosphate linkage**:
The unmodified linkage between two nucleotides.
_Avoid_: Normal bond, standard bond

## Modification

**Modification code**:
An IDT catalog identifier for a chemical group, written between slashes, e.g. `/5Biosg/`,
`/iSp18/`, `/3ThioMC3-D/`.
_Avoid_: Mod, tag, label

**Terminal modification**:
A modification code standing before the first or after the last nucleotide, so it belongs to
the 5' or 3' end of the chain.
_Avoid_: End mod

**Internal modification**:
A modification code standing between two nucleotides, so it belongs to the position it follows.
_Avoid_: Middle mod

**Thiol modifier**:
A pendant sulfhydryl on the end of a linker arm, `/5ThioMC6-D/`, `/3ThioMC3-D/` or
`/iThioMC6-D/`; a substituent for later coupling, not a backbone change.
_Avoid_: Thiolated linkage

**Sugar modification**:
A change to the ribose at a position: 2'-deoxy (DNA), ribose (`r`), 2'-O-methyl (`m`), LNA
(`+`), 2'-fluoro (`/i2FA/`, `/i2FC/`, `/i2FG/`, `/i2FU/`, `/i2FT/`) or 2'-MOE
(`/i2MOErA/`…`/i2MOErT/`, `/52MOEr…/`, `/32MOEr…/`).
_Avoid_: Base modification (that is the nucleobase, not the sugar)

**Degenerate base**:
A position holding a defined mixture of nucleobases, written with an IUB letter such as N, R or Y.
_Avoid_: Mixed base (reserve for the custom-mix form), ambiguity code

**Custom mix**:
A position written as parenthesised percentages in A C G T order, e.g. `(10203040)`.

## Drawing

**Schematic**:
The picture of the whole chain: tiles, connectors, tags, legend.
_Avoid_: Diagram, graphic, chart

**Tile**:
The coloured circle standing for one nucleotide.

**Sugar icon**:
The mini pentose ring drawn beside a structure base or below a linear tile, labelled `OH`, `OMe`,
`F`, `LNA` or `MOE`; absent for DNA, where the sugar is unmodified 2'-deoxyribose.
_Avoid_: Sugar glyph, corner letter

**Tag**:
The small labelled box showing a modification code verbatim, attached to the tile or the
connector it belongs to.
_Avoid_: Badge, chip, sticker

## Prediction

**Hairpin**:
A fold in which a single chain doubles back on itself: two antiparallel stem arms closed by a
loop. Distinct from a self-dimer, which involves two separate chains.

**Self-dimer**:
Two copies of the same oligonucleotide base-pairing with each other; IDT calls this a
homo-dimer. The tool shows the candidate with the lowest ΔG. Distinct from a hetero-dimer,
which involves two different sequences.
_Avoid_: Primer-dimer (may be hetero-dimer)

**Most likely homodimer**:
The IDT self-dimer candidate with the most negative ΔG under the stated Na⁺, Mg²⁺, oligo and
dNTP concentrations. The tool draws only this candidate.

**Oligo concentration**:
The strand concentration in µM used by the self-dimer calculation. It does not affect the
hairpin fold. IDT's qPCR default is 0.2 µM.

**dNTP concentration**:
The deoxynucleoside triphosphate concentration in mM used by the self-dimer calculation. It does
not affect the hairpin fold. IDT's qPCR default is 0.8 mM.

**Base-pair map**:
The list of which nucleotide positions pair with which, as predicted for one fold. It is not a
chemical drawing and does not show modification tags.

**Watson–Crick pair**:
A pair through the Watson–Crick edges: A–T, A–U or G–C. Drawn as a solid blue line.

**Wobble pair**:
A G·U (RNA) or G·T (DNA) pair through the Watson–Crick edges but with the bases displaced. Drawn
as a long-dashed teal line. It is a recognized non-Watson–Crick pair, not a generic mismatch.

**Noncanonical pair**:
A pair with a base identity other than Watson–Crick or G·U/G·T, such as G·A or A·C. Drawn as a
dotted amber line when IDT returns it. The drawing does not infer Hoogsteen or sugar-edge
geometry from the sequence alone.

**Hoogsteen edge / sugar edge**:
The other two base edges used in non-Watson–Crick pairing. Their geometry cannot be determined
from a two-dimensional base-pair map without structural data, so the tool does not use those
names for its line classes.

**Gibbs free energy (ΔG)**:
The predicted folding free energy in kcal/mol. More negative means a more favourable fold under
the stated conditions. The tool shows IDT's value; it never recomputes it.

**Melting temperature (Tm)**:
The temperature at which half the predicted structure is unfolded, in °C, as reported by IDT for
the chosen fold. The Analyzer property strip reports a separate oligo Tm from IDT's Analyze
endpoint; it is not the same as the chosen hairpin's Tm.

**Molecular weight**:
The oligo's calculated mass in g/mol, including its modifications, as reported by IDT.

**Extinction coefficient**:
The calculated A260 extinction coefficient in L/(mol·cm), including modifications, as reported
by IDT.

**OD260**:
One unit of absorbance at 260 nm in a 1 cm path. IDT reports how much oligo gives one OD260:
nmol/OD260 and µg/OD260.

**Lowest-ΔG fold**:
The structure with the most negative ΔG among all structures IDT returns, regardless of whether
it is a single hairpin or a more complex fold.

**Suboptimality**:
IDT's allowance, in percent, for returning structures above the minimum free energy. A larger
value returns more candidate folds; it does not change which returned fold is selected.

**Secondary-structure drawing**:
The view of the chosen fold: a stem-loop for one continuous hairpin, or a connected stem-loop
tree when the fold has several stems or a multiloop. A base-pair arc diagram is only a fallback.
2' sugar rings mark modified positions; `/Code/` modification tags remain in the linear view.

## The second meaning of "thiolated"

"Thiolated linkage" is ambiguous and is not used. Say **phosphorothioate linkage** for the
backbone sulfur, and **thiol modifier** for the pendant sulfhydryl.
