## Na3PS4 Reference Check

This example retains the supplied [POSCAR](POSCAR) and [reference lobsterin](reference.lobsterin)
unchanged. It contains 16 atoms: S 1–8, Na 9–14, and P 15–16. The two P centers form
finite PS4 units. One crosses cell boundaries, so cell translations matter.

From this folder, generate and render the configuration:

```bash
cat-prepare topology POSCAR --config seed.yaml -o config.yaml
cat-prepare lobsterin --config config.yaml -o generated.lobsterin
```

The [seed](seed.yaml) reproduces the reference's active requests, up to line ordering
and whitespace. It explicitly preserves the P–S and Na–S generator ranges.

| Requests | Count | Distances (Å) |
| --- | ---: | --- |
| P–S from generator | 8 | 2.059354 |
| Na–S from generator | 40 | 2.802646–3.411721 |
| S–S within PS4 units | 12 | 3.355268 or 3.378143 |
| S–S between PS4 units | 24 | 3.728470, 3.828392, or 3.872076 |

All 36 explicit atom-pair/cell requests match, including orbitalWise flags.
Pymatgen's independent periodic-neighbor enumeration also reproduces the complete
contact set. This is a geometry/input comparison; no LOBSTER calculation was run.

**Use `selection: all`, `cutoff: 4.2` for this reference.** With `nearest` and the
default 0.03 Å tolerance, only the 8 contacts at 3.728470 Å survive; the other 16
reference inter-component contacts are intentionally omitted. The next excluded
inter-component distance is 4.764656 Å, well beyond the chosen cutoff.

The reference includes basis entries for Sb and Se and four generators involving
those absent elements. The seed omits these entries because this POSCAR contains
only Na, P, and S. The active basis functions, energy window, projection/basis output
flags, and generator requests are preserved. Other compositions need their own seed.
