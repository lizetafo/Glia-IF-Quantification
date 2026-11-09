# Measurements

Column names start with `GFAP_` in the GFAP app and `IBA1_` in the IBA1 app; they are written as `MARKER_` below. Columns ending in `_um`, `_um2` or `_mm2` need a pixel size, which is read from `.im3` calibration. Without it they are `NA`.

## Processing
1. **Unmixing** of `.im3` cubes: reference library + NNLS, or library-free NMF. For composites, the channels are used directly.
2. **Background correction**: the marker image minus a Gaussian-blurred copy of itself (sigma = *Background sigma*), with negative values set to 0.
3. **Mask**: background-corrected marker > *threshold*, then objects smaller than *Remove objects smaller than* are dropped.
4. **Skeleton**: Zhang-Suen thinning of the mask, followed by the clean-up rules below.

## Area and intensity (primary outcomes)
| Column | Meaning |
|---|---|
| `MARKER_area_percent` | % of the field covered by the mask |
| `MARKER_positive_pixels`, `Total_pixels` | Mask pixels; all pixels |
| `MARKER_positive_area_um2`, `Full_FOV_area_um2` | The same, in µm² |
| `MARKER_mean_positive` | Mean marker signal inside the mask |
| `MARKER_mean_full_FOV` | Mean marker signal over the whole field |
| `MARKER_integrated_positive` | Summed marker signal inside the mask |
| `MARKER_signal_per_positive_um2_AU`, `..._full_FOV_um2_AU` | Summed signal per µm². Arbitrary units, because NMF signals are not calibrated intensities. |
| `DAPI_nuclei_rough` | Rough nucleus count (Otsu + watershed). Quality control only. |

## Skeleton (field-level descriptors)
| Column | Meaning |
|---|---|
| `MARKER_skeleton_length_px_rough` (`_um`, `_um_per_mm2`) | Weighted skeleton length: orthogonal step = 1, diagonal step = √2. A diagonal step is **not** counted where an orthogonal path already joins the two pixels. |
| `MARKER_internal_branchpoints_rough` (`_per_mm2`) | Skeleton pixels with ≥ 3 neighbours. Pixels **within 2 px of each other count as one branch point**. Pixels on the image border are excluded. |
| `MARKER_internal_endpoints_rough` (`_per_mm2`) | Skeleton pixels with exactly 1 neighbour, excluding the image border (where processes are cut off by the field edge) |
| `MARKER_skeleton_network_fragments_rough` | Number of separate skeleton pieces |
| `MARKER_area_to_skeleton_width_um_rough` | Mask area ÷ skeleton length: a rough mean process thickness |

**Skeleton clean-up (November 2026).** Zhang-Suen thinning leaves extra pixels at bends and "staircase" pixels on diagonal lines. These have 3 neighbours, so they were being counted as branch points, and they added extra length. They are now removed when they are not needed to keep the skeleton connected (Yokoi connectivity number = 1). Endpoints are never removed. See [validation](../validation/VALIDATION.md).

## Somas (IBA1 app only, rough)
Somas are found by **opening** the mask with a disc of radius *Soma opening radius*, which removes thin processes. Holes are then filled, touching somas are split by watershed, and objects outside *min–max soma area* are dropped.

| Column | Meaning |
|---|---|
| `IBA1_soma_objects_rough`, `IBA1_soma_density_per_mm2_rough` | Number of soma candidates, and per mm² |
| `IBA1_mean/median_soma_area_px_rough` (`_um2`) | Soma candidate size |
| `IBA1_mean_soma_intensity_AU_rough` | Mean IBA1 signal in soma candidates |
| `IBA1_soma_area_percent_of_FOV_rough`, `IBA1_process_area_percent_of_FOV_rough` | Field % covered by somas, and by the rest of the mask |
| `IBA1_soma_fraction_of_positive_area_percent_rough` | Somas as a % of the IBA1-positive area |
| `IBA1_..._per_soma_rough` | Field skeleton length, branch points and endpoints ÷ soma count. These are averages, not per-cell tracing. |

## Cautions
- All skeleton and soma values are **field-level**. Touching cells, 2-D projection of 3-D processes, the threshold, and GFAP-negative fine processes prevent true cell-by-cell reconstruction.
- Use **one locked threshold** for all images. Skeleton values depend strongly on the threshold.
- Results from before the November 2026 fix overcounted branch points 3–9× and length by 17–43%, so **do not pool them** with new results.
