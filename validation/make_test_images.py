"""
Synthetic ("phantom") glia images with exactly known answers.

Every structure is drawn from an ideal centerline, so the ground truth comes from
geometry, not from the analysis script being tested.

Outputs (in test_images/):
  1_topology_shapes.png    simple shapes: line, corner, T, Y, cross, loop, staircase
  2_microglia_phantom.png  IBA1-style cells (red) with soma + branched processes
  3_astrocyte_phantom.png  GFAP-style cells (red) around DAPI nuclei (blue), plus GFAP-negative nuclei
  4_microglia_realistic.png  same cells as 2 but blurred + noisy (robustness)
  5_crossings.png          four process crossings (X and + shapes)
  6_microglia_16bit_2ch.tif  cells of image 2 as a 16-bit, 2-channel, calibrated (0.5 um/px) TIFF
  ground_truth_*.csv
"""
import numpy as np, math, csv, os
from PIL import Image
from skimage.draw import line as skline, disk
from scipy.ndimage import binary_dilation, gaussian_filter
from scipy.spatial import ConvexHull

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "example_data", "synthetic")
os.makedirs(OUT, exist_ok=True)
ST = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)  # thickness: centerline dilated -> ~3 px wide


def draw_seg(m, p, q):
    rr, cc = skline(int(round(p[1])), int(round(p[0])), int(round(q[1])), int(round(q[0])))
    ok = (rr >= 0) & (rr < m.shape[0]) & (cc >= 0) & (cc < m.shape[1])
    m[rr[ok], cc[ok]] = True


def thick(center):
    return binary_dilation(center, ST)


def sholl_truth(segs, c, radii):
    """Exact number of crossings of each circle by the ideal centerline segments."""
    out = []
    for r in radii:
        n = 0
        for p, q in segs:
            dp, dq = math.dist(p, c), math.dist(q, c)
            lo, hi = min(dp, dq), max(dp, dq)
            # segments here are radially monotonic or start outside; count sign changes along segment
            ts = np.linspace(0, 1, 400)
            d = [math.dist((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])), c) - r for t in ts]
            n += sum(1 for a, b in zip(d[:-1], d[1:]) if (a < 0) != (b < 0))
        out.append(n)
    return out


def hull_area(mask):
    ys, xs = np.nonzero(mask)
    pts = np.concatenate([np.c_[xs, ys], np.c_[xs + 1, ys], np.c_[xs, ys + 1], np.c_[xs + 1, ys + 1]])
    return ConvexHull(pts).volume


def save_rgb(path, r, g=None, b=None):
    h, w = r.shape
    img = np.zeros((h, w, 3), np.uint8)
    img[..., 0] = r
    if g is not None: img[..., 1] = g
    if b is not None: img[..., 2] = b
    Image.fromarray(img).save(path)


# ---------------------------------------------------------------- 1. topology
H, W = 260, 700
center = np.zeros((H, W), bool)
shapes = []
def add(name, segs, ends, bps, branches):
    for p, q in segs: draw_seg(center, p, q)
    shapes.append((name, ends, bps, branches, sum(math.dist(p, q) for p, q in segs)))
add("straight line", [((20, 40), (90, 40))], 2, 0, 1)
add("right-angle corner", [((130, 20), (130, 70)), ((130, 70), (180, 70))], 2, 0, 1)
add("T junction", [((220, 30), (300, 30)), ((260, 30), (260, 90))], 3, 1, 3)
add("Y junction", [((380, 90), (380, 55)), ((380, 55), (350, 25)), ((380, 55), (410, 25))], 3, 1, 3)
add("cross", [((460, 50), (540, 50)), ((500, 15), (500, 90))], 4, 1, 4)
add("closed square loop", [((590, 20), (660, 20)), ((660, 20), (660, 90)), ((660, 90), (590, 90)), ((590, 90), (590, 20))], 0, 0, 1)
add("diagonal staircase", [((30, 150), (110, 190))], 2, 0, 1)
add("double branch", [((170, 220), (170, 140)), ((170, 190), (210, 160)), ((170, 165), (130, 135))], 4, 2, 5)
mask1 = thick(center)
save_rgb(f"{OUT}/1_topology_shapes.png", mask1.astype(np.uint8) * 200)
with open(f"{OUT}/ground_truth_1_topology.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["shape", "endpoints", "branch_points", "branches", "centerline_length_px"])
    for s in shapes: wr.writerow([s[0], s[1], s[2], s[3], round(s[4], 2)])
    wr.writerow(["TOTAL", sum(s[1] for s in shapes), sum(s[2] for s in shapes), sum(s[3] for s in shapes),
                 round(sum(s[4] for s in shapes), 2)])
print("topology totals: endpoints", sum(s[1] for s in shapes), "branch points", sum(s[2] for s in shapes))

# ---------------------------------------------------------------- 2. microglia
H, W = 520, 696
rng = np.random.default_rng(7)
SOMA_R = 7
cells = []
specs = [  # (cx, cy, [(angle_deg, length, secondary=(t, angle_offset, len) or None)])
    (120, 120, [(0, 60, None), (120, 55, (0.6, 40, 25)), (240, 60, None)]),
    (350, 110, [(30, 70, (0.5, -45, 30)), (100, 50, None), (190, 65, (0.55, 40, 28)), (280, 55, None)]),
    (580, 130, [(10, 50, None), (80, 45, None), (150, 55, None), (220, 50, None), (290, 45, None)]),
    (130, 380, [(45, 70, (0.4, 35, 30)), (135, 65, (0.5, -40, 25)), (225, 70, (0.45, 40, 30)), (315, 60, None)]),
    (370, 390, [(0, 45, None), (60, 50, None), (120, 45, None), (180, 50, None), (240, 45, None), (300, 50, None)]),
    (590, 390, [(90, 70, (0.5, 50, 30)), (270, 70, (0.5, -50, 30))]),
]
center = np.zeros((H, W), bool)
soma_mask = np.zeros((H, W), bool)
cell_masks = []
truth = []
for k, (cx, cy, procs) in enumerate(specs, 1):
    cm = np.zeros((H, W), bool)
    segs, primaries, ends, sec = [], 0, 0, 0
    for ang, L, s in procs:
        a = math.radians(ang)
        p0 = (cx + SOMA_R * math.cos(a), cy - SOMA_R * math.sin(a))
        p1 = (cx + (SOMA_R + L) * math.cos(a), cy - (SOMA_R + L) * math.sin(a))
        segs.append((p0, p1)); primaries += 1; ends += 1
        if s:
            t, off, L2 = s
            b0 = (p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1]))
            a2 = math.radians(ang + off)
            b1 = (b0[0] + L2 * math.cos(a2), b0[1] - L2 * math.sin(a2))
            segs.append((b0, b1)); ends += 1; sec += 1
    c1 = np.zeros((H, W), bool)
    for p, q in segs: draw_seg(c1, p, q)
    cm |= thick(c1)
    rr, cc = disk((cy, cx), SOMA_R + 0.5, shape=(H, W))
    sm = np.zeros((H, W), bool); sm[rr, cc] = True
    cm |= sm
    soma_mask |= sm
    cell_masks.append(cm)
    bp = (1 if primaries >= 3 else 0) + sec
    radii = list(range(5, 100, 5))
    sh = sholl_truth(segs, (cx + 0.5, cy + 0.5), radii)
    truth.append(dict(cell=k, soma_x=cx + 0.5, soma_y=cy + 0.5, primary_processes=primaries, endpoints=ends,
                      branch_points=bp, process_length_px=round(sum(math.dist(p, q) for p, q in segs), 2),
                      soma_area_px=int(sm.sum()), cell_area_px=int(cm.sum()), territory_area_px=round(hull_area(cm), 1),
                      sholl_max=max(sh), sholl=sh))
mask2 = np.any(cell_masks, axis=0)
red = np.where(mask2, 120, 0).astype(np.uint8)
red[soma_mask] = 220   # somas brighter, as for IBA1
save_rgb(f"{OUT}/2_microglia_phantom.png", red)
# realistic version: blur + noise
r2 = gaussian_filter(red.astype(float), 0.8) + rng.normal(0, 4, red.shape) + 6
save_rgb(f"{OUT}/4_microglia_realistic.png", np.clip(r2, 0, 255).astype(np.uint8))
with open(f"{OUT}/ground_truth_2_microglia.csv", "w", newline="") as f:
    cols = ["cell", "soma_x", "soma_y", "soma_area_px", "primary_processes", "endpoints", "branch_points",
            "process_length_px", "cell_area_px", "territory_area_px", "sholl_max"]
    wr = csv.writer(f); wr.writerow(cols)
    for t in truth: wr.writerow([t[c] for c in cols])
with open(f"{OUT}/ground_truth_2_microglia_sholl.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["cell", "radius_px", "intersections"])
    for t in truth:
        for r, n in zip(range(5, 100, 5), t["sholl"]): wr.writerow([t["cell"], r, n])
print("microglia: area % =", round(100 * mask2.mean(), 4), "cells", len(truth))

# ---------------------------------------------------------------- 3. astrocytes (GFAP + DAPI)
H, W = 520, 696
red = np.zeros((H, W), np.uint8); blue = np.zeros((H, W), float)
NUC_R = 6
astro = [(130, 130), (360, 140), (580, 120), (200, 380), (480, 370)]
negatives = [(60, 260), (300, 270), (640, 270), (360, 480)]
gmask = np.zeros((H, W), bool)
atruth = []
for k, (cx, cy) in enumerate(astro, 1):
    c1 = np.zeros((H, W), bool)
    n_proc = 5 + k % 3
    segs = []
    for j in range(n_proc):
        a = 2 * math.pi * j / n_proc + 0.3 * k
        p0 = (cx + (NUC_R + 2) * math.cos(a), cy - (NUC_R + 2) * math.sin(a))
        p1 = (cx + (NUC_R + 2 + 55) * math.cos(a), cy - (NUC_R + 2 + 55) * math.sin(a))
        segs.append((p0, p1))
    for p, q in segs: draw_seg(c1, p, q)
    # GFAP ring around the nucleus links the processes (cytoplasmic GFAP)
    rr, cc = disk((cy, cx), NUC_R + 3.5, shape=(H, W)); ring = np.zeros((H, W), bool); ring[rr, cc] = True
    rr, cc = disk((cy, cx), NUC_R + 1.5, shape=(H, W)); ring[rr, cc] = False
    cm = thick(c1) | ring
    gmask |= cm
    rr, cc = disk((cy, cx), NUC_R, shape=(H, W)); blue[rr, cc] = 180
    atruth.append(dict(cell=k, nucleus_x=cx + 0.5, nucleus_y=cy + 0.5, primary_processes=n_proc, endpoints=n_proc))
for (cx, cy) in negatives:
    rr, cc = disk((cy, cx), NUC_R, shape=(H, W)); blue[rr, cc] = 180
red[gmask] = 150
# GFAP fills the nucleus region weakly so positive cells have Cell: Red mean >= 4 (expanded cell region overlaps processes)
blue = np.clip(gaussian_filter(blue, 1.0), 0, 255).astype(np.uint8)
save_rgb(f"{OUT}/3_astrocyte_phantom.png", red, None, blue)
with open(f"{OUT}/ground_truth_3_astrocytes.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["measure", "value"])
    wr.writerow(["nuclei", len(astro) + len(negatives)])
    wr.writerow(["gfap_positive_cells", len(astro)])
    wr.writerow(["gfap_positive_cells_pct", round(100 * len(astro) / (len(astro) + len(negatives)), 2)])
    wr.writerow(["gfap_positive_area_pct", round(100 * gmask.mean(), 4)])
    wr.writerow(["gfap_positive_pixels", int(gmask.sum())])
    wr.writerow(["endpoints_total", sum(a["endpoints"] for a in atruth)])
    for a in atruth: wr.writerow([f"cell{a['cell']}_endpoints", a["endpoints"]])
print("astro: area % =", round(100 * gmask.mean(), 4), "pixels", int(gmask.sum()))

# ---------------------------------------------------------------- 5. crossings
m = np.zeros((200, 600), bool)
for p, q in [((20, 20), (120, 120)), ((20, 120), (120, 20)),       # diagonal X
             ((170, 70), (270, 70)), ((220, 20), (220, 120)),      # orthogonal +
             ((320, 30), (420, 110)), ((320, 110), (420, 40)),     # shallow X
             ((470, 20), (560, 180)), ((470, 180), (560, 30))]:    # steep X
    draw_seg(m, p, q)
save_rgb(f"{OUT}/5_crossings.png", thick(m).astype(np.uint8) * 200)
with open(f"{OUT}/ground_truth_5_crossings.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["measure", "value"]); wr.writerow(["endpoints", 16]); wr.writerow(["branch_points", 4])

# ---------------------------------------------------------------- 6. 16-bit multichannel, calibrated
import tifffile
r16 = np.array(Image.open(f"{OUT}/2_microglia_phantom.png"))[..., 0].astype(np.uint16) * 200
tifffile.imwrite(f"{OUT}/6_microglia_16bit_2ch.tif", np.stack([np.zeros_like(r16), r16]), imagej=True,
                 resolution=(2.0, 2.0), metadata={"axes": "CYX", "unit": "micron", "Labels": ["DAPI", "Iba1"]})
print("done")
