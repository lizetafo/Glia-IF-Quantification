# Validation

The skeleton measurements of both apps were tested **before and after** the November 2026 fix. Each test ran the functions inside the real `app.R` files under R 4.3 / EBImage 4.47. The GFAP and IBA1 apps share the same skeleton code and gave identical results.

Re-run everything from RStudio, in the repository folder:
```r
TEST_ORIGINAL <- TRUE   # also test the pre-fix versions
source("validation/test_morphology.R")
```
Results are written to `validation/morphology_results.csv`.

---

## 1. Images with known answers

Synthetic images drawn from ideal centerlines, in [`example_data/synthetic`](../example_data/synthetic). Each was thresholded at > 25.

| Image | True value | **After fix** | Before fix |
|---|---|---|---|
| Topology shapes: endpoints | 20 | **20** | 20 |
| Topology shapes: branch points | 5 | **5** | 23 |
| Topology shapes: length | 1134.3 px | **1121.5 (−1.1%)** | 1328.4 (+17%) |
| 6 microglia: endpoints / branch points | 32 / 13 | **32 / 13** | 32 / 40 |
| 6 microglia, blur + noise: endpoints / branch points | 32 / 13 | **32 / 13** | 32 / 69 |
| 4 crossings: endpoints / branch points | 16 / 4 | **16 / 4** | 16 / 35 |

After the fix, endpoints and branch points are **exact on all 8 measurements**. Before the fix, branch points were overcounted 3–9×.

Microglia skeleton length: 1774 px after the fix vs. 2537 px before (+43%). The validated QuPath implementation gives 1789 px on the same image.

## 2. Dense astrocyte-like networks vs. Fiji AnalyzeSkeleton

Four 696 × 520 fields with 25–55 overlapping branched cells ([`dense_fields/`](dense_fields)). Approximate true tip counts come from the drawing record ([`estimate_true_tips.py`](dense_fields/estimate_true_tips.py): drawn tips that don't overlap other processes and aren't on the border).

| Field | Endpoints: app / Fiji / ≈ true tips | Branch points: app / Fiji | Length (px): app / Fiji | Before fix: branch points / length |
|---|---|---|---|---|
| 1 | 204 / 242 / ≈178 | 211 / 244 | 6,641 / 7,221 | 410 / 9,644 |
| 2 | 214 / 268 / ≈178 | 284 / 331 | 7,843 / 8,494 | 468 / 11,034 |
| 3 | 253 / 318 / ≈213 | 363 / 429 | 10,110 / 10,856 | 641 / 14,695 |
| 4 | 324 / 423 / ≈248 | 508 / 600 | 12,868 / 13,885 | 812 / 18,388 |

- After the fix, the app gives **14–15% fewer** branch points, **16–23% fewer** endpoints and **7–8% less length** than Fiji. This is not caused by the border rule: only 0–2 endpoints per field lie on the border.
- The difference comes from the thinning algorithms. Fiji's thinning keeps more short twigs from noisy edges, while Zhang-Suen removes more of them. The app's endpoint counts are **closer to the true tip counts** than Fiji's, though both are above them.
- Before the fix, branch points were 35–68% **above** Fiji and length 30–35% above.

## 3. Real Iba1 images with hand annotations

Two images from the CellMorph dataset (Ling, 2026, CC BY 4.0; see [`example_data/cellmorph`](../example_data/cellmorph)). They were processed with the IBA1 app's own pipeline: background sigma 15, Otsu threshold, minimum object size 30 px, soma area 100–5000 px.

| Image | Soma overlap with hand outline (Dice, 0–1) | Hand-traced process points within 3 px of the skeleton | Median distance |
|---|---|---|---|
| CellMorph_001 | 0.57 (opening r = 3), 0.42 (r = 6) | 75% | 1.0 px |
| CellMorph_002 | 0.57 (r = 3), 0.53 (r = 6) | 53% | 2.2 px |

The skeleton follows hand-traced processes closely where the process is above the threshold. Faint stretches fall below the threshold and break the skeleton. Soma candidates only partly match the hand outlines, because bright neighbouring tissue merges into the soma. This supports the apps' description of soma values as **rough**, to be checked visually. Two images are too few for a general accuracy figure.

## 4. Full app run
Both apps were launched and driven through all tabs in a browser (upload → decompose → assign components → Otsu threshold → lock → batch → ZIP). On the GFAP + DAPI synthetic image, the GFAP app reported **5,152 positive pixels** and **9 nuclei**, both matching the true values exactly.

## Not tested
- `.im3` unmixing (NMF/NNLS): no `.im3` example file was available. Validate the component assignment against Nuance/inForm unmixed images.
- Pixel-size calibration from `.im3` metadata.
