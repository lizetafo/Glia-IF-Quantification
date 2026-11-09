import numpy as np
from skimage.draw import line
from skimage.filters import gaussian
from PIL import Image
rng = np.random.default_rng(42)
def astro_field(h=520, w=696, n_cells=25):
    img = np.zeros((h, w)); tips = []
    for _ in range(n_cells):
        cy, cx = rng.integers(30, h-30), rng.integers(30, w-30)
        for _ in range(rng.integers(4, 9)):
            ang = rng.uniform(0, 2*np.pi); L = rng.integers(20, 60)
            y1, x1 = int(cy+L*np.sin(ang)), int(cx+L*np.cos(ang))
            y1c, x1c = np.clip(y1, 0, h-1), np.clip(x1, 0, w-1)
            rr, cc = line(cy, cx, y1c, x1c); img[rr, cc] = 1; tips.append((y1c, x1c))
            for _ in range(rng.integers(0, 3)):
                t = rng.uniform(.4, .9); by, bx = int(cy+t*(y1-cy)), int(cx+t*(x1-cx))
                a2 = ang+rng.choice([-1, 1])*rng.uniform(.4, 1.2); L2 = rng.integers(8, 25)
                ey, ex = np.clip(int(by+L2*np.sin(a2)), 0, h-1), np.clip(int(bx+L2*np.cos(a2)), 0, w-1)
                rr, cc = line(np.clip(by, 0, h-1), np.clip(bx, 0, w-1), ey, ex); img[rr, cc] = 1; tips.append((ey, ex))
    centre = img.copy()
    img = gaussian(img, 1.5); img = img/img.max()*200 + rng.normal(0, 6, img.shape)
    return np.clip(img, 0, 255), centre, tips
for k in range(4):
    red, centre, tips = astro_field(n_cells=25+10*k)
    ref = np.array(Image.open(f'dense_field_{k+1}.png'))[..., 0]
    same = np.array_equal(red.astype(np.uint8), ref)
    # a true free tip: no other drawn centerline pixel within 4 px except along its own line (approx: count centre pixels in 9x9 window <= 5)
    free = 0
    for (y, x) in tips:
        win = centre[max(0, y-4):y+5, max(0, x-4):x+5]
        if win.sum() <= 5 and 0 < y < 519 and 0 < x < 695: free += 1
    print(f"dense_field_{k+1}: image reproduced={same}  drawn tips={len(tips)}  free (non-overlapping, not on border) tips={free}")
