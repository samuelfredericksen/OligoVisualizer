# OligoVisualizer V1

Paste an oligonucleotide in IDT's own notation and see a schematic of the molecule: one
horizontal 5' → 3' row, one coloured tile per nucleotide, the linkage drawn as the connector
between tiles, and every modification code drawn as a labelled tag.

The on-screen heading still reads "Oligo viewer"; this file is the project's documentation.

## Start it

Double-click **`Start Oligo Viewer.bat`**.

- The first run creates a private `.venv` folder and installs Flask from `requirements.txt`.
  That one step needs internet. Nothing is installed system-wide and no admin rights are used.
- The console prints `Oligo viewer running at http://127.0.0.1:8765/` and the browser opens.
  If 8765 is taken it tries upwards to 8785 and prints whichever port it got.
- Closing the console window stops the viewer. The page itself stores nothing, so nothing is
  lost by closing it.

## Read the drawing

|What you see|What it means|
|---|---|
|Coloured square with a letter|One nucleotide, in its base colour. Green A, blue C, amber G, vermilion T, pink U, light-blue I, grey for a degenerate base, pale grey `?` for a custom mix.|
|Small letter in the tile's top-right corner|The sugar. No letter = DNA (2'-deoxy). `r` = ribose, `m` = 2'-O-methyl, `+` = LNA. (`MOE` exists in the drawing code for 2'-MOE, but V1 cannot produce it — IDT writes 2'-MOE as a catalog modification, not a base prefix.)|
|Grey line between tiles|Ordinary phosphate linkage.|
|Dark red line with a circled **S**|Phosphorothioate linkage at that position. A short red stub after the last tile means the terminal linkage is thioated too.|
|Pale yellow box with a `/Code/`|A modification, printed exactly as you typed it. Above the first tile = 5' end, above the last tile = 3' end, between tiles = internal, attached to the nucleotide it follows.|
|Grey hatched tile with a letter|A degenerate base (N, R, Y, …).|
|Colour strip under `?`|A custom mix, split in A, C, G, T proportions.|
|`12 nt · 1 PS linkage · 2 modifications`|Length, number of thioated linkages, number of modification codes.|

A legend under the picture lists only the glyphs that sequence actually uses.

## What you can paste

```
/5Biosg/TTTTTT*GGGGGG/3Bio/           5' biotin, six thioated linkages, 3' biotin
/5ThioMC6-D/mA*+C*rG/3ThioMC3-D/      thiol modifier, 2'-OMe / LNA / ribose, two PS
ACGTN                                 plain DNA plus a degenerate last base
AC(10203040)T                         custom mix: 10% A, 20% C, 30% G, 40% T
```

Accepted: `A C G T U I` in any case; IUB codes `R Y S W K M B D H V N` in uppercase; `r`, `m`
and `+` prefixes; `*` after a nucleotide for a phosphorothioate; `/Code/` in any of IDT's
spellings; `(10203040)` or `(name:10203040)` for a custom mix. Pasted `5'-…-3'` quoting and
stray spaces or line breaks are ignored.

Two rules worth knowing:

- **Position comes from placement, not from the code's spelling.** `/ZEN/ACGT` draws ZEN as a 5'
  tag because it stands before the first base, even though the IDT catalog lists ZEN as an
  internal modification. Nothing is looked up: unknown and custom codes are drawn verbatim and
  never rejected.
- **`rT` and `mT` are rejected**, as IDT rejects them. You get a U tile plus a red note saying to
  write `rU`.

Bad input never blanks the picture: the valid part still draws and each problem appears below as
`position: message (token)`, in red for errors and amber for warnings. An unknown character
costs you that one character, not the sequence.

## Assumptions

- `A*C` and `A*` describe the same linkage; the `*` belongs to the nucleotide before it.
- Degenerate bases are hatched, custom mixes are split proportionally — a mixture with defined
  percentages is a different fact from an ambiguity code.
- Rounding is not involved anywhere: lengths and counts are exact.

## Limits

- About 500 nucleotides is the practical maximum; the drawing is one row and grows linearly, so
  beyond that it gets very wide and slow. There is no guard.
- This is a structure viewer only. Melting temperature, secondary structure, pKa and
  intramolecular interactions are not calculated here.

## Files

|File|What it is|
|---|---|
|`Start Oligo Viewer.bat`|The only thing you need to run. Sets up `.venv` once, then starts the viewer.|
|`app.py`|Local web server: serves the page and the `/api/render` endpoint (8765–8785).|
|`notation.py`|Parses the IDT string into a sequence model; collects problems instead of raising.|
|`render.py`|Draws the sequence model as SVG; owns the palette and geometry.|
|`static/index.html`|The page: textarea, debounced fetch, zoom, legend. No libraries, no build step.|
|`GLOSSARY.md`|Domain vocabulary, including why "thiolated" is not used.|
|`docs/adr/0001-server-rendered-svg.md`|Why the drawing happens in Python rather than the browser.|
|`tests/`|40 `unittest` tests: parser, renderer, HTTP contract.|

## Tests

```
py -m unittest tests.test_notation tests.test_render -v      # no venv needed
.venv\Scripts\python.exe -m unittest discover -v             # all 40, includes the server
```

## If something goes wrong

|Symptom|Meaning|
|---|---|
|Console says `Setup failed`|The one-time install could not download Flask. Check the network, then run the `.bat` again.|
|Page says `Cannot reach the viewer's server`|The console window was closed, or the file was opened directly instead of through the `.bat`. Start the `.bat` again and reload.|
|Console prints a port other than 8765|Something else already holds 8765. Use the URL the console printed.|
|No browser opened|Open the printed URL by hand; the viewer is running either way.|
