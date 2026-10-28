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
The mini pentose ring drawn below a tile, labelled `OH`, `OMe`, `F`, `LNA` or `MOE`; absent for
DNA, where the sugar is unmodified 2'-deoxyribose.
_Avoid_: Sugar glyph, corner letter

**Tag**:
The small labelled box showing a modification code verbatim, attached to the tile or the
connector it belongs to.
_Avoid_: Badge, chip, sticker

## The second meaning of "thiolated"

"Thiolated linkage" is ambiguous and is not used. Say **phosphorothioate linkage** for the
backbone sulfur, and **thiol modifier** for the pendant sulfhydryl.
