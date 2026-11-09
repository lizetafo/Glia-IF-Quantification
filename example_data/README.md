# Example data

## `cellmorph/`: real Iba1 microglia images
`CellMorph_001.png` and `CellMorph_002.png` (red = Iba1), each with a `.json` file containing a hand-drawn **soma outline** and the cell's **longest process path** (LabelMe format).

> Ling, H. (2026). *CellMorph dataset* (v1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.23585259 — licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Two of the 11 images are included here, unmodified apart from conversion from palette PNG to RGB PNG.

Open them in the **IBA1 app** with IBA1 = Red.

## `synthetic/`: images with known answers
Drawn from ideal centerlines, so the true endpoints, branch points and lengths are known exactly (`ground_truth_*.csv`). Used by `validation/test_morphology.R`.

| Image | Content | Use with |
|---|---|---|
| `1_topology_shapes.png` | line, corner, T, Y, cross, loop, staircase, double branch | either app (marker = Red) |
| `2_microglia_phantom.png` | 6 microglia-like cells | IBA1 app |
| `3_astrocyte_phantom.png` | 5 GFAP+ astrocytes (red) on DAPI nuclei (blue), plus 4 GFAP-negative nuclei | GFAP app (GFAP = Red, DAPI = Blue) |
| `4_microglia_realistic.png` | image 2 with blur and noise | IBA1 app |
| `5_crossings.png` | four process crossings | either app |

These images are released with this repository and may be reused freely.
