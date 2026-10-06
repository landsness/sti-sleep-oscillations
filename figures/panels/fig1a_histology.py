"""fig1a_histology.py -- Figure 1A: representative histology montage (2x2).

Columns STI+ / STI-; rows cresyl violet (cortical infarct) / GFAP (ipsilesional
thalamus). Raw TIFFs (no burned-in annotations) from data/histology/
(originals: DBSI Directory/histology/{cresyl_violet,IHC}/, also on Data@Becker).
Representative animals (locked): CV STI+ = M11 slice 8, CV STI- = M22 slice 13,
GFAP STI+ = M15 (WFCI_Ms15_2), GFAP STI- = M02 (Ms2_1).
All panels are mirrored so the (left-hemisphere) stroke is on the viewer's left.
Scale bars: 1 mm, calibrated per image from the Keyence BZ-X metadata embedded in each
TIFF (<Calibration>, nm/px; PlanApo 2x: 3.774 um/px at 1920x1440, 7.549 um/px for Ms2_1 at 960x720).

Output: results/panels/Fig1A_histology.png
"""
import io
import re
import struct
import sys
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont, ImageOps
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402

CV_POS = config.HISTOLOGY / "WFCI_Ms11_8.tif"
CV_NEG = config.HISTOLOGY / "WFCI_Ms22_13.tif"
GFAP_POS_SRC = config.HISTOLOGY / "WFCI_Ms15_2.tif"
GFAP_NEG_SRC = config.HISTOLOGY / "Ms2_1.tif"
OUT = config.PANELS / "Fig1A_histology.png"
FONT_BOLD = str(Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans-Bold.ttf")

def make_gfap_pos():
    im = ImageOps.mirror(Image.open(GFAP_POS_SRC).convert("RGB"))
    W,H = im.size
    arrows = [((26,44),(21,49.5)), ((16,58),(11.5,63.5))]  # mirrored, shortened
    fig = plt.figure(figsize=(W/100,H/100), dpi=100)
    ax = fig.add_axes([0,0,1,1]); ax.imshow(im); ax.axis("off")
    for (hx,hy),(tx,ty) in arrows:
        ax.annotate('', xy=(hx/100*W,hy/100*H), xytext=(tx/100*W,ty/100*H),
            arrowprops=dict(arrowstyle='-|>', color='white', lw=2.2, mutation_scale=24,
                shrinkA=0, shrinkB=1, path_effects=[pe.withStroke(linewidth=4, foreground='black')]))
    ax.set_xlim(0,W); ax.set_ylim(H,0)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=100); plt.close(fig); buf.seek(0)
    return Image.open(buf).convert("RGB")

def um_per_px(path):
    """Pixel size (um) from the Keyence BZ-X XML metadata block (<Calibration> is nm/px, stored as int64 bits of a double)."""
    s = Path(path).read_bytes().decode('latin1')
    m = re.search(r'<Calibration Type="System.Double">(-?\d+)</Calibration>', s)
    if not m:
        raise ValueError(f'no Keyence calibration in {path}')
    return struct.unpack('<d', struct.pack('<q', int(m.group(1))))[0] / 1000.0

TW, TH = 900, 675
FB = lambda s: ImageFont.truetype(FONT_BOLD, s)
def fit_scale(im):
    return max(TW/im.width, TH/im.height)
def fit(im):
    im = im.convert("RGB"); r = fit_scale(im)
    im = im.resize((int(im.width*r), int(im.height*r)))
    l=(im.width-TW)//2; t=(im.height-TH)//2
    return im.crop((l,t,l+TW,t+TH))
def sbar(t, white, px_per_mm):
    d=ImageDraw.Draw(t); col="white" if white else "black"
    d.rectangle([25,TH-45,25+round(px_per_mm),TH-33], fill=col)
    d.text((25,TH-82),"1 mm",fill=col,font=FB(28)); return t
def rot(txt,h):
    im=Image.new("RGB",(h,80),"white"); ImageDraw.Draw(im).text((10,20),txt,fill="black",font=FB(40))
    return im.rotate(90,expand=True)

def tile(src, img, white):
    """Crop/resize to the tile and draw a 1-mm bar using the source image's own pixel size."""
    px_per_mm = 1000.0 / um_per_px(src) * fit_scale(img)
    print(f'{Path(src).name}: {um_per_px(src):.3f} um/px -> 1 mm = {px_per_mm:.1f} tile px')
    return sbar(fit(img), white, px_per_mm)

t_cvpos = tile(CV_POS, ImageOps.mirror(Image.open(CV_POS)), False)
t_cvneg = tile(CV_NEG, ImageOps.mirror(Image.open(CV_NEG)), False)
t_gfpos = tile(GFAP_POS_SRC, make_gfap_pos(), True)
t_gfneg = tile(GFAP_NEG_SRC, ImageOps.mirror(Image.open(GFAP_NEG_SRC)), True)

LM,TMh,gap = 95,70,18
GW = LM+TW*2+gap*3; GH = TMh+TH*2+gap*3
M = Image.new("RGB",(GW,GH),"white"); d=ImageDraw.Draw(M)
d.text((LM+gap+TW//2-70,12),"STI+",fill="black",font=FB(46))
d.text((LM+gap*2+TW+TW//2-70,12),"STI−",fill="black",font=FB(46))
x1,x2=LM+gap,LM+gap*2+TW; y1,y2=TMh+gap,TMh+gap*2+TH
M.paste(t_cvpos,(x1,y1)); M.paste(t_cvneg,(x2,y1))
M.paste(t_gfpos,(x1,y2)); M.paste(t_gfneg,(x2,y2))
M.paste(rot("Cresyl violet",TH),(6,y1))
M.paste(rot("GFAP (thalamus)",TH),(6,y2))
d.text((10,4),"A",fill="black",font=FB(56))
OUT.parent.mkdir(parents=True, exist_ok=True)
M.save(OUT); print("saved", OUT, M.size)
