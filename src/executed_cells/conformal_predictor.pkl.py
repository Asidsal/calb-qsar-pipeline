# AUDIT RECORD for `conformal_predictor.pkl` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chemotype_assignments.csv
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
def P(name):
    """Locate a file by name in data/, results/ or models/ (intermediate files not shipped are noted in INDEX.md)."""
    for d in ("data", "results", "models"):
        p = ROOT / d / name
        if p.exists():
            return str(p)
    return str(ROOT / "results" / name)

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")

    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def set_frame(ax, style="open"):
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


def focal_palette(labels, focal, focal_color, other="muted", base_colors=None):
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    focal_set = {focal} if isinstance(focal, str) else set(focal)
    n = len(labels)
    if not focal_set & set(labels):
        raise ValueError(f"focal {focal!r} not found in labels")
    if base_colors is None:
        base_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#444444"])
    base_colors = [base_colors[i % len(base_colors)] for i in range(n)]
    if other == "grey":
        rest = ["#BCBCBC"] * n
    elif other == "ordinal":
        nf = max(1, n - len(focal_set))
        ramp = [mcolors.to_hex((v, v, v)) for v in
                ([0.55] if nf == 1 else [0.80 - 0.35 * i / (nf - 1) for i in range(nf)])]
        rest, k = [], 0
        for l in labels:
            rest.append(ramp[min(k, nf - 1)]); k += (l not in focal_set)
    else:
        def mute(c):
            r, g, b = mcolors.to_rgb(c)
            m = (r + g + b) / 3
            return mcolors.to_hex((0.3 * r + 0.7 * m, 0.3 * g + 0.7 * m, 0.3 * b + 0.7 * m))
        rest = [mute(c) for c in base_colors]
    return [focal_color if l in focal_set else rest[i] for i, l in enumerate(labels)]


def bar_with_points(ax, x, ymat, labels, colors, jitter=0.08, show_points=True,
                    errorbar=None, point_alpha=0.5, point_size=8):
    import numpy as np
    means = np.array([np.mean(y) for y in ymat], float)
    err = None
    if errorbar and not show_points:
        if errorbar == "sd":
            err = np.array([np.std(y, ddof=1) if np.asarray(y).size > 1 else 0 for y in ymat])
        elif errorbar == "ci95":
            from scipy.stats import t
            def _hw(y):
                n = np.asarray(y).size
                return t.ppf(0.975, n - 1) * np.std(y, ddof=1) / np.sqrt(n) if n > 1 else 0
            err = np.array([_hw(y) for y in ymat])
    ax.bar(x, means, color=colors, width=0.7, edgecolor="none",
           yerr=err, error_kw={"elinewidth": 0.8, "capsize": 0})
    if show_points:
        for xi, ys in zip(x, ymat):
            ys = np.asarray(ys)
            if ys.ndim and ys.size > 1:
                jit = (np.random.rand(ys.size) - 0.5) * 2 * jitter
                ax.scatter(np.full(ys.size, xi) + jit, ys, s=point_size, color="black",
                           alpha=point_alpha, zorder=3, linewidths=0)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    return ax


def strip_with_median(ax, groups, values, colors=None, jitter=0.12):
    import numpy as np
    labs = list(groups)
    if colors is None:
        colors = ["#444444"] * len(labs)
    for i, (ys, c) in enumerate(zip(values, colors)):
        ys = np.asarray(ys)
        jit = (np.random.rand(ys.size) - 0.5) * 2 * jitter
        ax.scatter(np.full(ys.size, i) + jit, ys, s=10, color=c, alpha=0.6, linewidths=0, zorder=2)
        m = np.median(ys)
        ax.plot([i - 0.22, i + 0.22], [m, m], color="black", lw=1.6, zorder=3)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs)
    return ax


def goodness_arrow(ax, text="higher = better", loc="upper left", axis="y", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams["legend.fontsize"]
    pos = {"upper left": (0.02, 0.98), "upper right": (0.98, 0.98),
           "lower left": (0.02, 0.02), "lower right": (0.98, 0.02)}[loc]
    ha = "left" if "left" in loc else "right"
    va = "top" if "upper" in loc else "bottom"
    arrow = "↑ " if axis == "y" else "→ "
    ax.text(pos[0], pos[1], arrow + text, transform=ax.transAxes,
            fontsize=fontsize, color=META_GREY, ha=ha, va=va)


def two_tier_label(name, meta):
    return f"{name}\n{meta}"


def end_of_line_labels(ax, xs, ys, labels, colors=None, dx=0.01, fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams["font.size"]
    if colors is None:
        colors = [None] * len(labels)
    span = ax.get_xlim()[1] - ax.get_xlim()[0]
    for x, y, lab, c in zip(xs, ys, labels, colors):
        ax.text(x[-1] + dx * span, y[-1], lab, color=c, va="center", ha="left", fontsize=fontsize)


import numpy as np, pandas as pd, json, time, warnings
warnings.filterwarnings('ignore')
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import matthews_corrcoef, roc_auc_score, balanced_accuracy_score, confusion_matrix
from xgboost import XGBClassifier

R = str(DATA_DIR)
SEED=42
Xtr_fp=np.load(R+r"\X_train_ecfp4.npy"); Xte_fp=np.load(R+r"\X_test_ecfp4.npy")
Xtr_d=pd.read_csv(R+r"\X_train_desc.csv"); Xte_d=pd.read_csv(R+r"\X_test_desc.csv")
ytr=np.load(R+r"\y_train.npy"); yte=np.load(R+r"\y_test.npy")
Xtr=np.hstack([Xtr_fp,Xtr_d.values]); Xte=np.hstack([Xte_fp,Xte_d.values])
sp=(ytr==0).sum()/ytr.sum()

def mk():
    sp=(ytr==0).sum()/ytr.sum()
    return XGBClassifier(n_estimators=500,learning_rate=0.05,max_depth=6,subsample=0.8,
        colsample_bytree=0.8,scale_pos_weight=sp,eval_metric='logloss',random_state=SEED,n_jobs=-1)

tr=pd.read_csv(R+r"\train.csv"); te=pd.read_csv(R+r"\test.csv")
chem=pd.read_csv(P("chemotype_assignments.csv"))[['std_smiles','chemotype']].drop_duplicates('std_smiles')
tr=tr.merge(chem,on='std_smiles',how='left'); te=te.merge(chem,on='std_smiles',how='left')

def build_cal(groups,a,y):
    return {(g,l):np.sort(a[(groups==g)&(y==l)]) for g in np.unique(groups) for l in (0,1)}

def pvals(P,groups,cal):
    pv=np.ones((len(P),2))
    for l in (0,1):
        a=1.0-P[:,l]
        for g in np.unique(groups):
            mk_=groups==g; ref=cal.get((g,l))
            if ref is None or len(ref)==0: continue
            pv[mk_,l]=(len(ref)-np.searchsorted(ref,a[mk_],side='left')+1)/(len(ref)+1)
    return pv

import numpy as np, pandas as pd, pickle, json
from sklearn.model_selection import train_test_split
idx=np.arange(len(ytr))
dp_prop,dp_cal=train_test_split(idx,test_size=0.25,stratify=ytr,random_state=SEED)
mdep=mk(); mdep.fit(Xtr[dp_prop],ytr[dp_prop])
Pc=mdep.predict_proba(Xtr[dp_cal]); yc=ytr[dp_cal]
ac=np.array([1-Pc[i,yc[i]] for i in range(len(yc))])
gc_raw=tr.chemotype.values[dp_cal]
cnt=pd.Series(gc_raw).value_counts(); keep_dep=set(cnt[cnt>=40].index)
gc=np.array([x if x in keep_dep else 'POOLED' for x in gc_raw])
cal_dep={(g,l):np.sort(ac[(gc==g)&(yc==l)]) for g in np.unique(gc) for l in (0,1)}
for l in (0,1): cal_dep[('POOLED',l)]=np.sort(ac[yc==l])
desc_cols=list(Xtr_d.columns)
pickle.dump(dict(model=mdep,cal=cal_dep,keep=keep_dep,desc_cols=desc_cols,
                 n_proper=len(dp_prop),n_cal=len(dp_cal),seed=SEED),open("conformal_predictor.pkl","wb"))