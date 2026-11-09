# Changelog

## November 2026: validated skeleton measurements
Applies to **both apps**.

**Fixed**
- **Branch points overcounted 3–9×.** Zhang-Suen thinning leaves extra pixels at bends and on diagonal "staircases". Each had 3 neighbours and was counted as a branch point. Such pixels are now removed when they are not needed for connectivity (`thin_redundant_pixels`). Junction pixels within 2 px of each other now count as one branch point (`count_components_within`).
- **Skeleton length overcounted by 17–43%.** At every bend, both the two short steps and the diagonal across them were added. A diagonal is now counted only where no orthogonal path joins the two pixels. Isolated single pixels no longer add length.

**Added**
- The GFAP app now finds Java automatically, as the IBA1 app does. Both apps also look in `C:\Program Files\Java\jre*`.
- A note in the Notes tab describing the skeleton rules.
- `setup/` scripts, user guide, measurement definitions, example data, and the validation suite.

**Unchanged**
Unmixing, background correction, thresholding, area and intensity measurements, DAPI counts, soma estimation, exports.

> Skeleton results produced before this version are not comparable with new results. Re-run the batch for all images.
