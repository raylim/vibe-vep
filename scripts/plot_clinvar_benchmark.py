"""Plot ClinVar benchmark results — 4-tool comparison."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Palette ───────────────────────────────────────────────────────────────────
BG      = "#1A1A2E"
PANEL   = "#22223A"
TEAL    = "#00C8AA"
AMBER   = "#FFC107"
CORAL   = "#FF6B6B"
LAVEND  = "#A78BFA"
GREY    = "#8888AA"
WHITE   = "#FFFFFF"
LGREY   = "#CCCCDD"

TOOLS   = ["vibe-vep", "snpEff", "VEP v115", "ANNOVAR"]
COLORS  = [TEAL, AMBER, CORAL, LAVEND]

fig = plt.figure(figsize=(18, 11), facecolor=BG)
fig.suptitle(
    "ClinVar Benchmark — 4-Tool Comparison\n"
    "232,008 pathogenic variants · independent ground truth · GENCODE v49 / Ensembl 115",
    color=WHITE, fontsize=15, fontweight="bold", y=0.98,
)

# ── Layout: 2×3 grid ─────────────────────────────────────────────────────────
gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.38,
                      left=0.06, right=0.97, top=0.90, bottom=0.07)

def style_ax(ax, title):
    ax.set_facecolor(PANEL)
    ax.spines[:].set_color("#333355")
    ax.tick_params(colors=LGREY, labelsize=9)
    ax.xaxis.label.set_color(LGREY)
    ax.yaxis.label.set_color(LGREY)
    ax.set_title(title, color=WHITE, fontsize=11, fontweight="bold", pad=8)

# ── 1. Overall HGVSp best vs any ──────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
style_ax(ax1, "Overall HGVSp Accuracy")

best = [87.4, 82.5, 79.8, 66.5]
any_ = [95.9, 95.9, 96.2, 93.4]
x    = np.arange(len(TOOLS))
w    = 0.35

bars1 = ax1.bar(x - w/2, best, w, color=COLORS, alpha=0.95, label="Best (primary tx)")
bars2 = ax1.bar(x + w/2, any_,  w, color=COLORS, alpha=0.45, label="Any transcript")

for bar, val in zip(bars1, best):
    ax1.text(bar.get_x() + bar.get_width()/2, val + 0.4,
             f"{val}%", ha="center", va="bottom", color=WHITE, fontsize=8, fontweight="bold")
for bar, val in zip(bars2, any_):
    ax1.text(bar.get_x() + bar.get_width()/2, val + 0.4,
             f"{val}%", ha="center", va="bottom", color=LGREY, fontsize=7.5)

ax1.set_xticks(x); ax1.set_xticklabels(TOOLS, fontsize=9)
ax1.set_ylim(55, 102); ax1.set_ylabel("% match", color=LGREY)
ax1.axhline(90, color=GREY, lw=0.6, ls="--", alpha=0.5)
solid = mpatches.Patch(color=GREY, alpha=0.95, label="Best (primary tx)")
light = mpatches.Patch(color=GREY, alpha=0.45, label="Any transcript")
ax1.legend(handles=[solid, light], fontsize=8, facecolor=PANEL,
           labelcolor=LGREY, edgecolor="#333355")

# ── 2. HGVSp by variant type ──────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
style_ax(ax2, "HGVSp by Variant Type (best)")

snv   = [90.2, 85.1, 81.2, None]
indel = [83.3, 78.6, 77.6, None]
x     = np.arange(len(TOOLS))

bars_s = ax2.bar(x - w/2, [v if v else 0 for v in snv],
                 w, color=COLORS, alpha=0.95, label="SNV")
bars_i = ax2.bar(x + w/2, [v if v else 0 for v in indel],
                 w, color=COLORS, alpha=0.45, label="Indel")

for bar, val in zip(bars_s, snv):
    if val:
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.4,
                 f"{val}%", ha="center", va="bottom", color=WHITE, fontsize=8, fontweight="bold")
    else:
        ax2.text(bar.get_x() + bar.get_width()/2, 1.5,
                 "N/A", ha="center", va="bottom", color=GREY, fontsize=7)
for bar, val in zip(bars_i, indel):
    if val:
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.4,
                 f"{val}%", ha="center", va="bottom", color=LGREY, fontsize=7.5)
    else:
        ax2.text(bar.get_x() + bar.get_width()/2, 1.5,
                 "N/A", ha="center", va="bottom", color=GREY, fontsize=7)

ax2.set_xticks(x); ax2.set_xticklabels(TOOLS, fontsize=9)
ax2.set_ylim(55, 100); ax2.set_ylabel("% match", color=LGREY)
solid = mpatches.Patch(color=GREY, alpha=0.95, label="SNV")
light = mpatches.Patch(color=GREY, alpha=0.45, label="Indel")
ax2.legend(handles=[solid, light], fontsize=8, facecolor=PANEL,
           labelcolor=LGREY, edgecolor="#333355")

# ── 3. Per-class best accuracy (grouped bar) ──────────────────────────────────
ax3 = fig.add_subplot(gs[0, 2])
style_ax(ax3, "HGVSp Best by Consequence Class")

classes   = ["missense", "frameshift", "stop\ngained", "inframe\ndel", "inframe\nins", "synonym."]
vibe_cls  = [91.3, 88.5, 83.1, 86.8, 83.3, 89.8]
snpeff_cls= [86.0, 84.1, 78.9, 82.1, 60.6,  0.0]
vep_cls   = [80.8, 82.3, 77.3, 78.3, 71.5,  0.0]
annov_cls = [70.6, 63.6, 63.0, 53.8,  0.0,  0.0]

xx  = np.arange(len(classes))
bw  = 0.20
off = [-1.5, -0.5, 0.5, 1.5]
all_vals = [vibe_cls, snpeff_cls, vep_cls, annov_cls]

for i, (vals, color) in enumerate(zip(all_vals, COLORS)):
    ax3.bar(xx + off[i]*bw, vals, bw, color=color, alpha=0.88, label=TOOLS[i])

ax3.set_xticks(xx); ax3.set_xticklabels(classes, fontsize=8)
ax3.set_ylim(0, 105); ax3.set_ylabel("% match", color=LGREY)
ax3.legend(fontsize=8, facecolor=PANEL, labelcolor=LGREY, edgecolor="#333355",
           loc="lower right")
ax3.axhline(90, color=GREY, lw=0.6, ls="--", alpha=0.5)

# ── 4. Failure taxonomy stacked bar ───────────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 0])
style_ax(ax4, "Failure Taxonomy")

transcript_err = [8.4, 13.4, 16.5, 30.6]
algo_fail      = [3.5,  3.9,  3.4,  4.3]
# not-annotated as % of total
not_ann_pct    = [1508/232008*100, 318/232008*100, 837/232008*100, 21/232008*100]

x = np.arange(len(TOOLS))
bw = 0.55

p1 = ax4.bar(x, transcript_err, bw, color=AMBER,  alpha=0.88, label="Transcript choice")
p2 = ax4.bar(x, algo_fail,      bw, color=CORAL,  alpha=0.88, label="Algorithmic fail",
             bottom=transcript_err)
p3 = ax4.bar(x, not_ann_pct,    bw, color=LAVEND, alpha=0.88, label="Not annotated",
             bottom=[t+a for t,a in zip(transcript_err, algo_fail)])

ax4.set_xticks(x); ax4.set_xticklabels(TOOLS, fontsize=9)
ax4.set_ylabel("% of variants", color=LGREY)
ax4.legend(fontsize=8, facecolor=PANEL, labelcolor=LGREY, edgecolor="#333355")

for xi, (t, a, n) in enumerate(zip(transcript_err, algo_fail, not_ann_pct)):
    total = t + a + n
    ax4.text(xi, total + 0.3, f"{total:.1f}%", ha="center", va="bottom",
             color=LGREY, fontsize=8)

# ── 5. Speed (log scale) ──────────────────────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 1])
style_ax(ax5, "Annotation Speed (log scale)")

speeds = [23226, 513, 174, 151]
x      = np.arange(len(TOOLS))
bars   = ax5.bar(x, speeds, 0.55, color=COLORS, alpha=0.88)
ax5.set_yscale("log")
ax5.set_xticks(x); ax5.set_xticklabels(TOOLS, fontsize=9)
ax5.set_ylabel("variants / second", color=LGREY)
ax5.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
    lambda v, _: f"{int(v):,}"))

for bar, val, tool in zip(bars, speeds, TOOLS):
    ax5.text(bar.get_x() + bar.get_width()/2,
             val * 1.25,
             f"{val:,} v/s", ha="center", va="bottom",
             color=WHITE if tool == "vibe-vep" else LGREY,
             fontsize=8, fontweight="bold" if tool == "vibe-vep" else "normal")

# speedup annotations
for i, (s, tool) in enumerate(zip(speeds[1:], TOOLS[1:]), 1):
    ratio = speeds[0] / s
    ax5.text(i, s * 0.35, f"{ratio:.0f}× slower",
             ha="center", va="top", color=GREY, fontsize=7.5, style="italic")

# ── 6. Consequence class match ────────────────────────────────────────────────
ax6 = fig.add_subplot(gs[1, 2])
style_ax(ax6, "Consequence Class Match")

conseq = [97.7, 97.0, 96.9, 97.8]
x      = np.arange(len(TOOLS))
bars   = ax6.bar(x, conseq, 0.55, color=COLORS, alpha=0.88)
ax6.set_xticks(x); ax6.set_xticklabels(TOOLS, fontsize=9)
ax6.set_ylim(95.5, 98.5); ax6.set_ylabel("% match", color=LGREY)

for bar, val in zip(bars, conseq):
    ax6.text(bar.get_x() + bar.get_width()/2, val + 0.04,
             f"{val}%", ha="center", va="bottom", color=WHITE, fontsize=9, fontweight="bold")

ax6.axhline(97, color=GREY, lw=0.6, ls="--", alpha=0.5)

# ── Save ──────────────────────────────────────────────────────────────────────
out = "testdata/clinvar/benchmark_plot.png"
fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=BG)
print(f"Saved: {out}")
