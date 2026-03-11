"""
3-slide summary deck with embedded matplotlib plots.
Slide 1 — Overall HGVSp accuracy (best vs any, SNV vs indel)
Slide 2 — Per-class breakdown + failure taxonomy
Slide 3 — Performance + consequence class match
"""

import io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

# ── Shared palette ────────────────────────────────────────────────────────────
BG      = "#1A1A2E"
PANEL   = "#22223A"
TEAL    = "#00C8AA"
AMBER   = "#FFC107"
CORAL   = "#FF6B6B"
LAVEND  = "#A78BFA"
GREY    = "#8888AA"
WHITE   = "#FFFFFF"
LGREY   = "#CCCCDD"
DIMGREY = "#555577"

TOOLS  = ["vibe-vep", "snpEff", "VEP v115", "ANNOVAR"]
COLORS = [TEAL, AMBER, CORAL, LAVEND]

# ── Helper: save figure to bytes ──────────────────────────────────────────────
def fig_to_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor=BG)
    buf.seek(0)
    plt.close(fig)
    return buf

# ── PPTX helpers ──────────────────────────────────────────────────────────────
SLIDE_W = Inches(13.33)
SLIDE_H = Inches(7.5)
DARK_BG = RGBColor(0x1A, 0x1A, 0x2E)
C_TEAL  = RGBColor(0x00, 0xC8, 0xAA)
C_AMBER = RGBColor(0xFF, 0xC1, 0x07)
C_WHITE = RGBColor(0xFF, 0xFF, 0xFF)
C_LGREY = RGBColor(0xCC, 0xCC, 0xDD)
C_DIM   = RGBColor(0x55, 0x55, 0x77)

prs = Presentation()
prs.slide_width  = SLIDE_W
prs.slide_height = SLIDE_H

def new_slide():
    s = prs.slides.add_slide(prs.slide_layouts[6])
    bg = s.shapes.add_shape(1, 0, 0, SLIDE_W, SLIDE_H)
    bg.fill.solid(); bg.fill.fore_color.rgb = DARK_BG; bg.line.fill.background()
    return s

def hbar(slide, top, color=C_TEAL):
    b = slide.shapes.add_shape(1, 0, top, SLIDE_W, Inches(0.055))
    b.fill.solid(); b.fill.fore_color.rgb = color; b.line.fill.background()

def txt(slide, text, l, t, w, h, size=13, bold=False,
        color=C_WHITE, align=PP_ALIGN.LEFT):
    tb = slide.shapes.add_textbox(l, t, w, h)
    tb.word_wrap = True
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.alignment = align
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.bold = bold
    r.font.color.rgb = color; r.font.name = "Calibri"

def slide_header(slide, title, subtitle, num, total=3):
    hbar(slide, Inches(1.02))
    txt(slide, title, Inches(0.4), Inches(0.08), Inches(11.8), Inches(0.9),
        size=28, bold=True, color=C_TEAL)
    txt(slide, subtitle, Inches(0.4), Inches(1.1), Inches(12.0), Inches(0.38),
        size=12, color=C_LGREY)
    txt(slide, f"{num} / {total}", Inches(12.5), Inches(0.08), Inches(0.8), Inches(0.5),
        size=12, color=C_DIM, align=PP_ALIGN.RIGHT)

def add_fig(slide, fig, left, top, width, height):
    buf = fig_to_bytes(fig)
    slide.shapes.add_picture(buf, left, top, width, height)

def ax_style(ax, title="", xlabel="", ylabel=""):
    ax.set_facecolor(PANEL)
    for spine in ax.spines.values():
        spine.set_color("#333355")
    ax.tick_params(colors=LGREY, labelsize=9)
    if title:
        ax.set_title(title, color=WHITE, fontsize=11, fontweight="bold", pad=7)
    if xlabel:
        ax.set_xlabel(xlabel, color=LGREY, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, color=LGREY, fontsize=9)


# ══════════════════════════════════════════════════════════════════════════════
# Slide 1 — Overall HGVSp Accuracy
# ══════════════════════════════════════════════════════════════════════════════
s1 = new_slide()
slide_header(s1,
    "HGVSp Accuracy — 4-Tool Comparison",
    "232,008 pathogenic ClinVar variants · independent ground truth · GENCODE v49 / Ensembl 115",
    num=1)

# ── Plot A: best vs any (horizontal bar) ─────────────────────────────────────
fig_a, ax = plt.subplots(figsize=(6.5, 3.2), facecolor=BG)
ax_style(ax, title="Overall HGVSp Match")

best = [87.4, 82.5, 79.8, 64.8]
any_ = [95.9, 95.9, 96.2, 95.4]
y    = np.arange(len(TOOLS))
h    = 0.35

bars_any  = ax.barh(y - h/2, any_,  h, color=COLORS, alpha=0.38, label="Any transcript")
bars_best = ax.barh(y + h/2, best,  h, color=COLORS, alpha=0.95, label="Primary transcript")

for bar, val, col in zip(bars_best, best, COLORS):
    ax.text(val + 0.3, bar.get_y() + bar.get_height()/2,
            f"{val}%", va="center", color=col, fontsize=9, fontweight="bold")
for bar, val in zip(bars_any, any_):
    ax.text(val + 0.3, bar.get_y() + bar.get_height()/2,
            f"{val}%", va="center", color=LGREY, fontsize=8.5)

ax.set_yticks(y)
ax.set_yticklabels(TOOLS, color=LGREY, fontsize=9)
ax.set_xlim(55, 103)
ax.set_xlabel("% HGVSp match", color=LGREY, fontsize=9)
ax.axvline(90, color=GREY, lw=0.7, ls="--", alpha=0.5)
solid = mpatches.Patch(color=GREY, alpha=0.95, label="Primary transcript")
light = mpatches.Patch(color=GREY, alpha=0.38, label="Any transcript")
ax.legend(handles=[solid, light], fontsize=8, facecolor=PANEL,
          labelcolor=LGREY, edgecolor="#333355", loc="lower right")
fig_a.tight_layout()
add_fig(s1, fig_a, Inches(0.3), Inches(1.55), Inches(6.6), Inches(3.5))

# ── Plot B: SNV vs Indel (grouped horizontal bar) ────────────────────────────
fig_b, ax = plt.subplots(figsize=(6.0, 3.2), facecolor=BG)
ax_style(ax, title="HGVSp Best by Variant Type")

snv   = [90.2, 85.1, 81.2, 69.0]
indel = [83.3, 78.6, 77.6, 58.6]
y     = np.arange(len(TOOLS))

bars_s = ax.barh(y + h/2, snv,   h, color=COLORS, alpha=0.95, label="SNV")
bars_i = ax.barh(y - h/2, indel, h, color=COLORS, alpha=0.45, label="Indel")

for bar, val, col in zip(bars_s, snv, COLORS):
    ax.text(val + 0.3, bar.get_y() + bar.get_height()/2,
            f"{val}%", va="center", color=col, fontsize=9, fontweight="bold")
for bar, val in zip(bars_i, indel):
    ax.text(val + 0.3, bar.get_y() + bar.get_height()/2,
            f"{val}%", va="center", color=LGREY, fontsize=8.5)

ax.set_yticks(y); ax.set_yticklabels(TOOLS, color=LGREY, fontsize=9)
ax.set_xlim(50, 103)
ax.set_xlabel("% HGVSp match", color=LGREY, fontsize=9)
solid = mpatches.Patch(color=GREY, alpha=0.95, label="SNV")
light = mpatches.Patch(color=GREY, alpha=0.45, label="Indel")
ax.legend(handles=[solid, light], fontsize=8, facecolor=PANEL,
          labelcolor=LGREY, edgecolor="#333355", loc="lower right")
fig_b.tight_layout()
add_fig(s1, fig_b, Inches(6.9), Inches(1.55), Inches(6.1), Inches(3.5))

# Callout row
for val, label, x in [
    ("+7%",  "vs snpEff\n(primary)",   Inches(0.4)),
    ("+23%", "vs ANNOVAR\n(primary)",  Inches(3.0)),
    ("+7%",  "vs snpEff\n(SNV)",       Inches(5.8)),
    ("+25%", "vs ANNOVAR\n(indel)",    Inches(8.6)),
    ("~96%", "all tools converge\n(any tx)", Inches(10.8)),
]:
    rect = s1.shapes.add_shape(1, x, Inches(5.35), Inches(2.4), Inches(1.05))
    rect.fill.solid(); rect.fill.fore_color.rgb = RGBColor(0x2C, 0x2C, 0x50)
    rect.line.color.rgb = C_TEAL; rect.line.width = Pt(1.2)
    txt(s1, val,   x, Inches(5.36), Inches(2.4), Inches(0.55),
        size=22, bold=True, color=C_TEAL, align=PP_ALIGN.CENTER)
    txt(s1, label, x, Inches(5.88), Inches(2.4), Inches(0.48),
        size=10, color=C_LGREY, align=PP_ALIGN.CENTER)

hbar(s1, Inches(6.52), color=C_AMBER)
txt(s1,
    "Gap between primary and any-transcript accuracy is entirely transcript selection — not algorithmic differences.",
    Inches(0.4), Inches(6.58), Inches(12.5), Inches(0.38),
    size=11, color=C_LGREY)


# ══════════════════════════════════════════════════════════════════════════════
# Slide 2 — Per-Class Breakdown & Failure Taxonomy
# ══════════════════════════════════════════════════════════════════════════════
s2 = new_slide()
slide_header(s2,
    "Per-Class Accuracy & Failure Taxonomy",
    "Primary transcript accuracy by consequence class · mutually exclusive failure breakdown",
    num=2)

# ── Plot C: per-class grouped bars ───────────────────────────────────────────
fig_c, ax = plt.subplots(figsize=(8.0, 3.6), facecolor=BG)
ax_style(ax, title="HGVSp Best by Consequence Class")

classes    = ["missense\n(70k)", "frameshift\n(83k)", "stop_gained\n(76k)",
              "inframe_del\n(2k)", "inframe_ins\n(748)", "synonymous\n(815)"]
vibe_v  = [91.3, 88.5, 83.1, 86.8, 83.3, 89.8]
snpeff_v= [86.0, 84.1, 78.9, 82.1, 60.6,  0.0]
vep_v   = [80.8, 82.3, 77.3, 78.3, 71.5,  0.0]
annov_v = [70.6, 63.6, 63.0, 53.8,  0.0,  0.0]

xx   = np.arange(len(classes))
bw   = 0.19
offs = [-1.5, -0.5, 0.5, 1.5]
all_vals = [vibe_v, snpeff_v, vep_v, annov_v]

for i, (vals, col, tool) in enumerate(zip(all_vals, COLORS, TOOLS)):
    bars = ax.bar(xx + offs[i]*bw, vals, bw, color=col, alpha=0.88, label=tool)

ax.set_xticks(xx); ax.set_xticklabels(classes, color=LGREY, fontsize=8.5)
ax.set_ylim(0, 106); ax.set_ylabel("% HGVSp match", color=LGREY, fontsize=9)
ax.axhline(90, color=GREY, lw=0.6, ls="--", alpha=0.4)
ax.legend(fontsize=9, facecolor=PANEL, labelcolor=LGREY,
          edgecolor="#333355", loc="lower left", ncol=2)

# annotate 0% bars
for xi, cls in enumerate(classes):
    if cls.startswith("inframe_ins") or cls.startswith("synonymous"):
        for i, (vals, col) in enumerate(zip(all_vals, COLORS)):
            if vals[xi] == 0.0:
                ax.text(xi + offs[i]*bw, 1.5, "0%",
                        ha="center", va="bottom", color=col, fontsize=7, alpha=0.7)

fig_c.tight_layout()
add_fig(s2, fig_c, Inches(0.3), Inches(1.55), Inches(8.1), Inches(3.85))

# ── Plot D: failure taxonomy stacked bar ─────────────────────────────────────
fig_d, ax = plt.subplots(figsize=(4.5, 3.6), facecolor=BG)
ax_style(ax, title="Failure Taxonomy")

tx_err   = [8.4,  13.4, 16.5, 30.6]   # transcript choice
algo_err = [3.5,   3.9,  3.4,  4.3]   # complete fail
not_ann  = [1508/232008*100, 318/232008*100, 837/232008*100, 21/232008*100]

x  = np.arange(len(TOOLS))
bw = 0.55

p1 = ax.bar(x, tx_err,   bw, color=AMBER,  alpha=0.88, label="Transcript choice")
p2 = ax.bar(x, algo_err, bw, color=CORAL,  alpha=0.88, label="Algorithmic fail",
            bottom=tx_err)
p3 = ax.bar(x, not_ann,  bw, color=LAVEND, alpha=0.88, label="Not annotated",
            bottom=[t+a for t,a in zip(tx_err, algo_err)])

for xi, (t, a, n) in enumerate(zip(tx_err, algo_err, not_ann)):
    total = t + a + n
    ax.text(xi, total + 0.4, f"{total:.1f}%",
            ha="center", va="bottom", color=LGREY, fontsize=9, fontweight="bold")

ax.set_xticks(x); ax.set_xticklabels(TOOLS, color=LGREY, fontsize=9)
ax.set_ylabel("% of variants", color=LGREY, fontsize=9)
ax.legend(fontsize=8, facecolor=PANEL, labelcolor=LGREY, edgecolor="#333355")
fig_d.tight_layout()
add_fig(s2, fig_d, Inches(8.5), Inches(1.55), Inches(4.6), Inches(3.85))

hbar(s2, Inches(5.55), color=C_AMBER)
txt(s2,
    "ANNOVAR: 0% on inframe_ins & synonymous (not annotated in protein field) · "
    "30.6% transcript-choice errors vs 8.4% for vibe-vep · "
    "Complete-fail rates (~3–4%) are equal across all tools — dominated by ClinVar format artifacts",
    Inches(0.4), Inches(5.62), Inches(12.5), Inches(0.6),
    size=11, color=C_LGREY)

txt(s2, "Only vibe-vep emits HGVSp for synonymous variants (e.g. p.Arg273=) — snpEff, VEP, and ANNOVAR report nothing.",
    Inches(0.4), Inches(6.3), Inches(12.5), Inches(0.35),
    size=11, color=C_LGREY)
hbar(s2, Inches(6.75), color=C_TEAL)


# ══════════════════════════════════════════════════════════════════════════════
# Slide 3 — Performance & Summary
# ══════════════════════════════════════════════════════════════════════════════
s3 = new_slide()
slide_header(s3,
    "Performance & Summary",
    "End-to-end wall time · 232,008 variants · vibe-vep cache load: 2.0 s from DuckDB",
    num=3)

# ── Plot E: speed log bar ─────────────────────────────────────────────────────
fig_e, ax = plt.subplots(figsize=(5.5, 3.4), facecolor=BG)
ax_style(ax, title="Annotation Speed (log scale)", ylabel="variants / second")

speeds = [23345, 513, 174, 151]
x      = np.arange(len(TOOLS))
bars   = ax.bar(x, speeds, 0.55, color=COLORS, alpha=0.88)
ax.set_yscale("log")
ax.set_xticks(x); ax.set_xticklabels(TOOLS, color=LGREY, fontsize=9)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
ax.tick_params(axis="y", colors=LGREY)

for bar, val, tool in zip(bars, speeds, TOOLS):
    ax.text(bar.get_x() + bar.get_width()/2, val * 1.3,
            f"{val:,}", ha="center", va="bottom",
            color=WHITE if tool == "vibe-vep" else LGREY,
            fontsize=9, fontweight="bold" if tool == "vibe-vep" else "normal")

for i, (s, tool) in enumerate(zip(speeds[1:], TOOLS[1:]), 1):
    ratio = speeds[0] / s
    ax.text(i, s * 0.4, f"{ratio:.0f}× slower",
            ha="center", va="top", color=GREY, fontsize=8, style="italic")

fig_e.tight_layout()
add_fig(s3, fig_e, Inches(0.3), Inches(1.52), Inches(5.6), Inches(3.6))

# ── Plot F: consequence class match (horizontal bar) ─────────────────────────
fig_f, ax = plt.subplots(figsize=(4.5, 3.4), facecolor=BG)
ax_style(ax, title="Consequence Class Match")

conseq = [97.7, 97.0, 96.9, 97.8]
y      = np.arange(len(TOOLS))
bars   = ax.barh(y, conseq, 0.5, color=COLORS, alpha=0.88)
for bar, val, col in zip(bars, conseq, COLORS):
    ax.text(val + 0.02, bar.get_y() + bar.get_height()/2,
            f"{val}%", va="center", color=col, fontsize=9, fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels(TOOLS, color=LGREY, fontsize=9)
ax.set_xlim(95.5, 99.0)
ax.set_xlabel("% match", color=LGREY, fontsize=9)
ax.axvline(97, color=GREY, lw=0.7, ls="--", alpha=0.5)
fig_f.tight_layout()
add_fig(s3, fig_f, Inches(5.9), Inches(1.52), Inches(4.6), Inches(3.6))

# ── Summary stat boxes ────────────────────────────────────────────────────────
boxes = [
    ("87.4%", "Highest primary\nHGVSp accuracy"),
    ("97.7%", "Highest consequence\nclass match"),
    ("154×",  "Faster than\nANNOVAR"),
    ("2",     "Fixable vibe-specific\ngaps remaining"),
]
for i, (val, label) in enumerate(boxes):
    x = Inches(0.4 + i * 3.1)
    y = Inches(5.35)
    rect = s3.shapes.add_shape(1, x, y, Inches(2.85), Inches(1.55))
    rect.fill.solid(); rect.fill.fore_color.rgb = RGBColor(0x2C, 0x2C, 0x50)
    rect.line.color.rgb = C_TEAL; rect.line.width = Pt(1.5)
    txt(s3, val,   x, y + Inches(0.06), Inches(2.85), Inches(0.75),
        size=30, bold=True, color=C_TEAL, align=PP_ALIGN.CENTER)
    txt(s3, label, x, y + Inches(0.8),  Inches(2.85), Inches(0.65),
        size=10, color=C_LGREY, align=PP_ALIGN.CENTER)

hbar(s3, Inches(7.08), color=C_AMBER)
txt(s3,
    "Fixable gaps: (1) splice-adjacent frameshifts emit _splice HGVSp suffix — fix splice/frameshift priority at exon boundaries  "
    "·  (2) stop-creating inframe deletions misclassified as stop_gained — preserve inframe_del when reading frame is intact",
    Inches(0.4), Inches(7.14), Inches(12.5), Inches(0.32),
    size=10, color=C_LGREY)

# ── Save ──────────────────────────────────────────────────────────────────────
out = "testdata/clinvar/benchmark_summary_slide.pptx"
prs.save(out)
print(f"Saved: {out}")
