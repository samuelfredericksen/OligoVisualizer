# OligoVisualizer 2.0

Paste an oligonucleotide in IDT's own notation and see a schematic of the molecule: one
horizontal 5' → 3' row, one coloured circle per nucleotide, a mini sugar ring below each
modified position, the linkage drawn as the connector between circles, and every modification
code drawn as a labelled tag.

The **Secondary structure** tab asks IDT's OligoAnalyzer for hairpin foldings, picks the one with
the lowest Gibbs free energy (ΔG), and draws its base-pair map as native SVG. The linear view is
still the default and still works without an IDT account.

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
|Coloured circle with a letter|One nucleotide, in its base colour. Green A, blue C, amber G, vermilion T, pink U, light-blue I, grey for a degenerate base, pale grey `?` for a custom mix.|
|Mini pentose ring below a circle|The sugar at that position. No ring = DNA (2'-deoxy). The ring's label reads `OH` = ribose, `OMe` = 2'-O-methyl, `F` = 2'-fluoro, `LNA` = locked nucleic acid (blue bridge), `MOE` = 2'-O-methoxyethyl.|
|Grey line between circles|Ordinary phosphate linkage.|
|Dark red line with a circled **S**|Phosphorothioate linkage at that position. A short red stub after the last circle means the terminal linkage is thioated too.|
|Pale yellow box with a `/Code/`|A modification, printed exactly as you typed it. Above the first circle = 5' end, above the last circle = 3' end, between circles = internal, attached to the nucleotide it follows.|
|Grey hatched circle with a letter|A degenerate base (N, R, Y, …).|
|Colour strip inside `?`|A custom mix, split in A, C, G, T proportions.|
|`12 nt · 1 PS linkage · 2 modifications`|Length, number of thioated linkages, number of modification codes.|

A legend under the picture lists only the glyphs that sequence actually uses.

## Predict a secondary structure

1. Click **Sign in to IDT** and enter your IDT account in the local form. The password is sent to
   this local viewer and forwarded to IDT over HTTPS; it is not stored or logged. Only the IDT
   session and anti-forgery tokens stay in the viewer process's memory.
2. Click **Predict structure**. The tool calls the real IDT OligoAnalyzer hairpin endpoint with
   the same default conditions the Analyzer uses for hairpin folding: DNA mode, 25 °C, 50 mM
   Na⁺, 3 mM Mg²⁺, 50% suboptimality, 20 maximum foldings, linear sequence. DNA is IDT's
   default mode even when the sequence contains U.
3. The view switches to the fold with the lowest **ΔG** in kcal/mol (not just the first
   single-hairpin fold). The summary shows ΔG, Tm and which structure was chosen; the detail line
   shows ΔH, ΔS and the conditions.
4. A fold that is one continuous stem-loop is drawn as circles on two arms with a loop. A fold
   with multiple stems is drawn as connected stem-loops: separate hairpins sit side by side,
   nested stems stack, and multiloops branch above their parent stem. A labelled arc diagram is
   only a last-resort fallback. The structure drawing shows bases only; modification codes remain
   visible in the linear view.

Modified oligos are sent to IDT exactly as written. The tool does not strip `*`, `+`, `r`, `m`,
`/Code/` tags or custom mixes; the Analyzer decides whether it accepts a particular chemistry.
If it rejects one, the IDT error is shown and the linear view remains available.

Predictions are cached per sequence for the running viewer. Switching tabs reuses the cached
result, and clicking **Predict structure** again reuses it without a second IDT call. **Sign out**
clears the IDT session; predictions then ask for sign-in again. The password field is cleared
immediately after submission.

## What you can paste

```
/5Biosg/TTTTTT*GGGGGG/3Bio/           5' biotin, six thioated linkages, 3' biotin
/5ThioMC6-D/mA*+C*rG/3ThioMC3-D/      thiol modifier, 2'-OMe / LNA / ribose, two PS
/i2FA/*/i2FC/*/i2FG/*/i2FU/           2'-fluoro A/C/G/U, three thioated linkages
/52MOErA/*/i2MOErC/*/32MOErT/         2'-MOE bases, internal and terminal, two PS
ACGTN                                 plain DNA plus a degenerate last base
AC(10203040)T                         custom mix: 10% A, 20% C, 30% G, 40% T
```

Accepted: `A C G T U I` in any case; IUB codes `R Y S W K M B D H V N` in uppercase; `r`, `m`
and `+` prefixes; IDT 2'-F base codes (`/i2FA/`, `/i2FC/`, `/i2FG/`, `/i2FU/`, `/i2FT/`) and
2'-MOE base codes (`/i2MOErA/`…`/i2MOErT/`, `/52MOErA/`…`/52MOErT/`,
`/32MOErA/`…`/32MOErT/`); `*` after a nucleotide for a phosphorothioate; `/Code/` in any of
IDT's spellings; `(10203040)` or `(name:10203040)` for a custom mix. Pasted `5'-…-3'` quoting and
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
- IDT 2'-F and 2'-MOE base codes are nucleotides, not modification codes: they draw a circle
  with the code's base letter and a sugar icon, and do not count as modifications in the summary.
- Degenerate bases are hatched, custom mixes are split proportionally — a mixture with defined
  percentages is a different fact from an ambiguity code.
- Lengths and counts are exact. ΔG, Tm, ΔH and ΔS are IDT's values, formatted for display only;
  the tool never recomputes them.

## Limits

- About 500 nucleotides is the practical maximum for the linear drawing; it is one row and grows
  linearly, so beyond that it gets very wide and slow. There is no guard.
- Secondary structure requires internet access and a signed-in IDT account. The prediction uses
  IDT's defaults listed above; they are displayed but not editable in this version.
- ΔG and Tm are IDT's values, not recalculated here. Complex folds are drawn as connected
  stem-loops, not as a full two-dimensional chemical diagram with sugar and modification details.
- If your IDT account presents an MFA/CAPTCHA challenge, this local form cannot complete it. Use
  IDT in a normal browser for that account.

## Files

|File|What it is|
|---|---|
|`Start Oligo Viewer.bat`|The only thing you need to run. Sets up `.venv` once, then starts the viewer.|
|`app.py`|Local web server: page, `/api/render`, `/api/structure`, and `/api/idt/*` (8765–8785).|
|`idt.py`|IDT sign-in, hairpin request, CT base-pair map retrieval, ΔG selection.|
|`notation.py`|Parses the IDT string into a sequence model; collects problems instead of raising.|
|`render.py`|Draws the linear schematic and secondary-structure SVG; owns palette and geometry.|
|`static/index.html`|The page: textarea, IDT sign-in, tabs, debounced fetch, zoom, legend. No build step.|
|`dev_run.py`|Optional development launcher; reads gitignored `.env` for auto-login. Not used by the `.bat`.|
|`GLOSSARY.md`|Domain vocabulary, including why "thiolated" is not used.|
|`docs/adr/0001-server-rendered-svg.md`|Why the drawing happens in Python rather than the browser.|
|`tests/`|107 `unittest` tests: parser, renderers, IDT client, HTTP contract.|

## Development auto-login (optional)

For repeated local testing only, create a gitignored `.env` beside `app.py`:

```
IDT_USERNAME=your-idt-username
IDT_PASSWORD=your-idt-password
```

Then run `.venv\Scripts\python.exe dev_run.py`. It signs in to IDT once and serves at
`http://127.0.0.1:8767/`. The normal `Start Oligo Viewer.bat` ignores `.env`; nothing is installed
system-wide. `.gitignore` is not encryption, so remove the file when you no longer need it.

## Tests

```
py -m unittest tests.test_notation tests.test_render tests.test_idt -v   # no network, no venv
.venv\Scripts\python.exe -m unittest discover -v                        # all 107, server too
```

No test contacts idtdna.com; the IDT client is exercised with fake HTTP methods and fake
responses.

## If something goes wrong

|Symptom|Meaning|
|---|---|
|Console says `Setup failed`|The one-time install could not download Flask. Check the network, then run the `.bat` again.|
|Page says `Cannot reach the viewer's server`|The console window was closed, or the file was opened directly instead of through the `.bat`. Start the `.bat` again and reload.|
|Console prints a port other than 8765|Something else already holds 8765. Use the URL the console printed.|
|No browser opened|Open the printed URL by hand; the viewer is running either way.|
|Secondary structure says `Sign in to IDT`|Sign in with the local form; the session lives only in the running viewer process.|
|IDT sign-in or prediction gives an HTTP/auth error|IDT may be down, the login page may have changed, or the session expired. Sign in again; the linear view is unaffected.|
|A sequence has no predicted structure|IDT found no base pairs. This is a warning, not a drawing failure.|
|MFA/CAPTCHA appears on IDT's site|This local form cannot complete it. Use IDT in a normal browser for that account.|
