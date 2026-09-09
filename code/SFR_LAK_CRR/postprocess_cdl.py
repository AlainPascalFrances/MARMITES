"""
================================================================================
Post-processing for the CdL groundwater-flow model (cdl_gwf_model_opusv1.py)
================================================================================
Reads the MODFLOW 6 output already on disk (no need to re-run the model) and
produces, in WORKSPACE/postproc/:

  1. Head TIME SERIES at the obs_points (cdl_gwf.obs.head.csv)   -> figure + tidy CSV
  2. Mean-head MAPS per layer                                    -> figure + CSV
  2b. Mean-DEPTH MAPS per layer (land surface − head)            -> figure + CSV
  3. WATER BUDGET by compartment (surface / unsaturated / aquifer):
       - overall mean rate (m3/d) from the listing budget        -> figure + CSV
       - yearly volumes (m3/yr)                                   -> figure + CSV
       - unsaturated-zone internal budget (UZF .bud)             -> figure + CSV
       - surface-water internal budget (SFR .bud)                -> figure + CSV
       - per-layer storage change (GWF .cbc)                     -> figure + CSV

The spin-up (steady-state period 0 + the 11 transient months that replay year 1)
is EXCLUDED from every average/map so the statistics reflect the real period.

Run in the activated env, e.g.:
  & "C:\\miniconda3\\Scripts\\conda.exe" run -p C:\\miniconda3\\envs\\flopy ^
      --no-capture-output python -u postprocess_cdl.py
================================================================================
"""

import re
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects

import flopy

# --- CONFIG (keep in sync with cdl_gwf_model_opusv1.py) ----------------------
WORKSPACE   = Path(r"E:\00code_ws\DRYAD\CdL_model")
MODEL_NAME  = "cdl_gwf"
SIM_START   = pd.Timestamp(1981, 1, 1)
SPINUP_NPER = 12          # SS period 0 + 11 transient spin-up months -> first 12 SPs
GPKG        = r"E:/zzCloud/OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia/DRYAD/GIS/dryad_modelo_NbS.gpkg"
OBS_LAYER   = "obs_points_cdl"
P_CSV       = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\p_month_198101_202605.csv"
ET_CSV      = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\et0_month_198101_202605.csv"

PIEZO_XLSX  = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\piezos_qualitative_month_198101_202605.xlsx"

# Match the model run's start stamp (written by the main script) so preproc\<stamp>\
# and postproc\<stamp>\ correspond to the same run; fall back to now if absent.
_stampf = WORKSPACE / "last_run_stamp.txt"
RUN_STAMP = _stampf.read_text().strip() if _stampf.exists() else datetime.now().strftime("%Y%m%d%H%M")
OUT = WORKSPACE / "postproc" / RUN_STAMP
OUT.mkdir(exist_ok=True, parents=True)

HDS  = WORKSPACE / f"{MODEL_NAME}.hds"
CBC  = WORKSPACE / f"{MODEL_NAME}.cbc"
LST  = WORKSPACE / f"{MODEL_NAME}.lst"
UZFB = WORKSPACE / f"{MODEL_NAME}.uzf.bud"
SFRB = WORKSPACE / f"{MODEL_NAME}.sfr.bud"
OBSCSV = WORKSPACE / f"{MODEL_NAME}.obs.head.csv"
GRB  = WORKSPACE / f"{MODEL_NAME}.disv.grb"
ROUTED_SHP = WORKSPACE / "streams_cdl_routed.shp"


# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------
def load_grid():
    """Return (modelgrid, top[ncpl], nlay) from the binary grid file (fallback: pkl)."""
    try:
        grb = flopy.mf6.utils.MfGrdFile(str(GRB))
        mg = grb.modelgrid
        return mg, np.asarray(mg.top, dtype=float).ravel(), mg.nlay
    except Exception as e:
        print(f"   (grb load failed: {e!r}; falling back to voronoi_grid.pkl, no top)")
        import pickle
        with open(WORKSPACE / "voronoi_grid.pkl", "rb") as f:
            gridprops_vg, _ = pickle.load(f)
        # nlay inferred later from heads; build a 1-layer grid just for geometry
        from flopy.discretization import VertexGrid
        mg = VertexGrid(**gridprops_vg, nlay=1)
        return mg, None, None


def real_kstpkper(binfile):
    """(kstp,kper) saved AFTER the spin-up (kper >= SPINUP_NPER)."""
    return [kk for kk in binfile.get_kstpkper() if kk[1] >= SPINUP_NPER]


def real_dates(n):
    """Monthly calendar dates for the n real stress periods (post spin-up)."""
    return pd.date_range(SIM_START, periods=n, freq="MS")


def read_forcing(csv, value_col="MONTHLY_TOTAL"):
    """Monthly totals (mm) as a Series indexed by month-start date (same parse as the model)."""
    df = pd.read_csv(csv, sep=r"\s+|,|;", engine="python")
    datec = [c for c in df.columns if "date" in c.lower()][0]
    if value_col not in df.columns:                       # fall back to last column
        value_col = df.columns[-1]
    df[datec] = pd.to_datetime(df[datec], format="%m/%d/%Y", errors="coerce")
    s = (df[[datec, value_col]].rename(columns={datec: "date", value_col: "val"})
         .dropna().set_index("date").sort_index())
    return s["val"]


def compartment_of(term):
    t = term.upper()
    if "STO" in t:
        return "aquifer (storage)"
    if "UZF" in t:
        return "unsaturated (UZF->GW)"
    if "SFR" in t:
        return "surface (stream)"
    if "DRN" in t:
        return "surface (drain)"
    return "other"


COMP_COLORS = {
    "aquifer (storage)": "tab:brown",
    "unsaturated (UZF->GW)": "tab:green",
    "surface (stream)": "tab:blue",
    "surface (drain)": "tab:cyan",
    "other": "0.6",
}

# Brown palette for the head curves — ties to the aquifer/storage color (tab:brown)
# and keeps the curves distinct from the blue rainfall / orange ET0 bars.
HEAD_COLORS = ["#5c3211", "#8c564b", "#b3743a", "#cf9d62"]  # dark brown -> tan (#8c564b = tab:brown)


def read_piezo_obs():
    """Observed heads per piezometer from the Excel (one sheet per point, p0..p6).
    Returns {POINTNAME: DataFrame[date, head]} keyed to match the obs points (P0..).
    These are the PEST calibration targets (qualitative for now)."""
    out = {}
    try:
        xl = pd.ExcelFile(PIEZO_XLSX)
    except Exception as e:
        print(f"   (piezo Excel not read: {e!r})")
        return out
    for sh in xl.sheet_names:
        m = re.match(r"\s*(p\d+)", sh, re.IGNORECASE)
        name = m.group(1).upper() if m else sh
        d = pd.read_excel(PIEZO_XLSX, sheet_name=sh)
        datec = [c for c in d.columns if "date" in c.lower()]
        headc = [c for c in d.columns if "head" in c.lower()]
        if not datec or not headc:
            continue
        s = (d[[datec[0], headc[0]]].rename(columns={datec[0]: "date", headc[0]: "head"}))
        s["date"] = pd.to_datetime(s["date"], errors="coerce")
        out[name] = s.dropna().sort_values("date")
    return out


# =============================================================================
# 1. HEAD TIME SERIES AT OBS POINTS
# =============================================================================
def obs_timeseries():
    print(">> [1] Obs-point head time series …")
    if not OBSCSV.exists():
        print(f"   !! {OBSCSV.name} not found — was the OBS package written/run? Skipping.")
        return
    df = pd.read_csv(OBSCSV)
    tcol = df.columns[0]
    t = df[tcol].to_numpy(dtype=float)

    # drop the spin-up portion of the record (model time < end of spin-up)
    spin = pd.date_range(SIM_START, periods=SPINUP_NPER + 1, freq="MS")
    spinup_days = (spin[-1] - spin[0]).days
    keep = t >= spinup_days - 1e-6
    dates = SIM_START + pd.to_timedelta(t[keep] - spinup_days, unit="D")

    # group columns by point name (strip _L<k>)
    pts = {}
    for c in df.columns[1:]:
        m = re.match(r"(.+)_L(\d+)$", c, re.IGNORECASE)
        name, lay = (m.group(1), int(m.group(2))) if m else (c, None)
        pts.setdefault(name, []).append((lay, c))

    # tidy long CSV
    tidy = []
    for name, cols in pts.items():
        for lay, c in cols:
            for d, v in zip(dates, df.loc[keep, c].to_numpy()):
                tidy.append((name, lay, d, v))
    pd.DataFrame(tidy, columns=["point", "layer", "date", "head_m"]) \
        .to_csv(OUT / "obs_head_timeseries.csv", index=False)

    # figure: one panel per point, a line per layer + land-surface reference
    mg, top, _ = load_grid()
    try:
        xc = np.array(mg.xcellcenters).ravel()
        yc = np.array(mg.ycellcenters).ravel()
    except Exception:
        xc = yc = None
    pt_node = {}
    if xc is not None:
        import geopandas as gpd
        from scipy.spatial import cKDTree
        gobs = gpd.read_file(GPKG, layer=OBS_LAYER).to_crs(3763)
        tree = cKDTree(np.column_stack([xc, yc]))
        ncol = "Name" if "Name" in gobs.columns else gobs.columns[0]
        for _, r in gobs.iterrows():
            if r.geometry is not None and not r.geometry.is_empty:
                _, idx = tree.query([r.geometry.x, r.geometry.y])
                pt_node[str(r[ncol]).strip().replace(" ", "_")] = int(idx)

    # monthly rainfall & ET0 over the plotted span -> bars hanging from a top axis
    months = pd.date_range(dates.min().to_period("M").to_timestamp(),
                           dates.max().to_period("M").to_timestamp(), freq="MS")
    rain = et = None
    try:
        _r = read_forcing(P_CSV).reindex(months).to_numpy(dtype=float)
        _e = read_forcing(ET_CSV).reindex(months).to_numpy(dtype=float)
        pe_max = float(np.nanmax(np.concatenate([_r, _e])))
        if not np.isfinite(pe_max):
            raise ValueError("no finite rainfall/ET in the plotted span")
        rain, et = _r, _e
    except Exception as e:
        print(f"   (rainfall/ET bars skipped: {e!r})")

    obs_meas = read_piezo_obs()    # observed heads per point (overlay + later PEST targets)
    names = list(pts.keys())
    n = len(names)
    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 2.9 * nrows),
                             squeeze=False, sharex=True)
    for ax, name in zip(axes.ravel(), names):
        for j, (lay, c) in enumerate(sorted(pts[name])):
            col = HEAD_COLORS[((lay - 1) if lay else j) % len(HEAD_COLORS)]
            ax.plot(dates, df.loc[keep, c].to_numpy(), lw=1.3, color=col, zorder=5,
                    label=f"L{lay}" if lay else c)
        if top is not None and name in pt_node:
            ax.axhline(top[pt_node[name]], color="0.5", ls="--", lw=0.8,
                       label="land surface", zorder=4)
        if name in obs_meas and len(obs_meas[name]):       # observed heads (PEST targets)
            om = obs_meas[name]
            _m = (om["date"] >= dates.min()) & (om["date"] <= dates.max())
            ax.scatter(om.loc[_m, "date"], om.loc[_m, "head"], s=26, marker="o",
                       facecolor="k", edgecolor="w", linewidths=0.4, zorder=8,
                       label="observed")
        ax.set_ylabel("head (m)")
        ax2 = None
        if rain is not None:
            ax2 = ax.twinx()
            w = 12.0                                  # bar width (days); two per month
            mid = months + pd.Timedelta(days=14)
            ax2.bar(mid - pd.Timedelta(days=w / 2), rain, width=w,
                    color="tab:blue", alpha=0.75, label="Rainfall")
            ax2.bar(mid + pd.Timedelta(days=w / 2), et, width=w,
                    color="tab:orange", alpha=0.75, label="ET$_0$")
            ax2.set_ylim(pe_max * 3.2, 0.0)           # 0 at top -> bars descend
            ax2.set_ylabel("P, ET$_0$ (mm/mo)", fontsize=8)
            ax2.tick_params(labelsize=7)
            ax.set_zorder(ax2.get_zorder() + 1)       # head lines in front of bars
            ax.patch.set_visible(False)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels() if ax2 is not None else ([], [])
        ax.legend(h1 + h2, l1 + l2, fontsize=7, ncol=3, loc="lower left")
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.3)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    fig.suptitle("CdL — simulated head at observation points "
                 "(rainfall & ET$_0$ on top; post spin-up)", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "obs_head_timeseries.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote obs_head_timeseries.png / .csv  ({n} points)")


# =============================================================================
# 2. MEAN-HEAD & MEAN-DEPTH MAPS PER LAYER
# =============================================================================
def compute_mean_head():
    """Mean head per layer (nlay, ncpl) over the post-spin-up periods, + kstpkper used."""
    hds = flopy.utils.HeadFile(str(HDS))
    kk = real_kstpkper(hds) or hds.get_kstpkper()
    arr = np.array([np.squeeze(hds.get_data(kstpkper=k)) for k in kk])  # (nt,nlay,ncpl)
    if arr.ndim == 2:                       # single layer -> (nt,ncpl)
        arr = arr[:, None, :]
    arr = np.where(np.abs(arr) > 1e29, np.nan, arr)
    return np.nanmean(arr, axis=0), kk      # (nlay,ncpl), list


def load_overlays():
    """(obs_points gdf, routed-streams gdf) for map overlays; None where unavailable."""
    gobs = routed = None
    try:
        import geopandas as gpd
        gobs = gpd.read_file(GPKG, layer=OBS_LAYER).to_crs(3763)
    except Exception:
        pass
    try:
        import geopandas as gpd
        if ROUTED_SHP.exists():
            routed = gpd.read_file(ROUTED_SHP)
    except Exception:
        pass
    return gobs, routed


def _plot_layer_maps(arrays, titles, cbar_label, cmap, suptitle, outfile,
                     mg, gobs, routed):
    """Shared per-layer map figure (one panel per array) used by head & depth maps."""
    n = len(arrays)
    ncols = min(n, 3)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows), squeeze=False)
    axf = axes.ravel()
    for i, (a, t) in enumerate(zip(arrays, titles)):
        ax = axf[i]
        pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax, layer=0)
        ca = pmv.plot_array(a, cmap=cmap)
        pmv.plot_grid(lw=0.1, color="0.85")
        if routed is not None:
            routed.plot(ax=ax, color="white", lw=0.8, zorder=3)
        if gobs is not None:
            _ncol = "Name" if "Name" in gobs.columns else gobs.columns[0]
            ax.scatter(gobs.geometry.x, gobs.geometry.y, marker="s", s=40,
                       c="magenta", edgecolors="k", linewidths=0.6, zorder=5)
            for _, _r in gobs.iterrows():
                if _r.geometry is None or _r.geometry.is_empty:
                    continue
                ax.annotate(str(_r[_ncol]), (_r.geometry.x, _r.geometry.y),
                            textcoords="offset points", xytext=(5, 4),
                            fontsize=8, fontweight="bold", color="k", zorder=6,
                            path_effects=[patheffects.withStroke(linewidth=2,
                                                                 foreground="white")])
        fig.colorbar(ca, ax=ax, shrink=0.7, label=cbar_label)
        ax.set_title(t)
        ax.set_aspect("equal")
    for ax in axf[n:]:
        ax.set_visible(False)
    fig.suptitle(suptitle, y=1.0)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def head_maps():
    print(">> [2] Mean-head maps per layer …")
    if not HDS.exists():
        print(f"   !! {HDS.name} not found. Skipping.")
        return
    mean_head, kk = compute_mean_head()
    nlay = mean_head.shape[0]
    mg, top, _ = load_grid()
    gobs, routed = load_overlays()
    pd.DataFrame({"layer": np.arange(1, nlay + 1),
                  "mean_head_m": np.nanmean(mean_head, axis=1)}) \
        .to_csv(OUT / "mean_head_per_layer.csv", index=False)
    _plot_layer_maps(
        [mean_head[l] for l in range(nlay)],
        [f"Mean head — layer {l + 1}" for l in range(nlay)],
        "head (m)", "viridis",
        f"CdL — mean head over {len(kk)} post-spin-up periods",
        OUT / "mean_head_maps.png", mg, gobs, routed)
    print(f"   wrote mean_head_maps.png / mean_head_per_layer.csv  ({nlay} layers)")


def depth_maps():
    """Per-layer mean DEPTH = land surface (topography) − head, mirroring head_maps."""
    print(">> [2b] Mean-depth maps per layer (land surface − head) …")
    if not HDS.exists():
        print(f"   !! {HDS.name} not found. Skipping.")
        return
    mg, top, _ = load_grid()
    if top is None:
        print(f"   !! no top elevation ({GRB.name} missing) — cannot compute depth. Skipping.")
        return
    mean_head, kk = compute_mean_head()
    nlay = mean_head.shape[0]
    gobs, routed = load_overlays()
    depth = np.array([top - mean_head[l] for l in range(nlay)])   # + = head below ground
    pd.DataFrame({"layer": np.arange(1, nlay + 1),
                  "mean_depth_m": np.nanmean(depth, axis=1)}) \
        .to_csv(OUT / "mean_depth_per_layer.csv", index=False)
    _plot_layer_maps(
        [depth[l] for l in range(nlay)],
        [f"Mean depth (surface − head) — layer {l + 1}" for l in range(nlay)],
        "depth below ground (m)", "RdYlBu_r",
        f"CdL — mean depth to head over {len(kk)} post-spin-up periods\n"
        f"(+ = head below land surface,  − = above)",
        OUT / "mean_depth_maps.png", mg, gobs, routed)
    print(f"   wrote mean_depth_maps.png / mean_depth_per_layer.csv  ({nlay} layers)")


# =============================================================================
# 3a. COMPARTMENT BUDGET FROM THE LISTING (overall mean + yearly)
# =============================================================================
def list_budget():
    print(">> [3a] Compartment water budget (listing) …")
    if not LST.exists():
        print(f"   !! {LST.name} not found. Skipping.")
        return
    mflist = flopy.utils.Mf6ListBudget(str(LST))
    df_flux, _ = mflist.get_dataframes(start_datetime=str(SIM_START.date()), diff=False)
    if df_flux is None or len(df_flux) == 0:
        print("   !! no budget parsed. Skipping.")
        return
    flux = df_flux.iloc[SPINUP_NPER:].copy()           # drop SS + spin-up
    n_real = len(flux)
    flux.index = real_dates(n_real)

    drop = {"TOTAL", "PERCENT_DISCREPANCY", "IN-OUT"}
    bases = sorted({c[:-3] for c in flux.columns
                    if c.endswith("_IN") and c[:-3] not in drop})
    net = pd.DataFrame(
        {b: flux.get(b + "_IN", 0.0) - flux.get(b + "_OUT", 0.0) for b in bases},
        index=flux.index)                              # + = source INTO aquifer

    # ---- overall mean rate (m3/d) ----
    summ = pd.DataFrame({
        "term": bases,
        "compartment": [compartment_of(b) for b in bases],
        "mean_rate_m3d": [net[b].mean() for b in bases],
    }).sort_values(["compartment", "term"])
    summ.to_csv(OUT / "budget_compartment_mean.csv", index=False)

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = [COMP_COLORS.get(c, "0.6") for c in summ["compartment"]]
    ax.barh(summ["term"], summ["mean_rate_m3d"], color=colors, edgecolor="k", lw=0.4)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean rate (m³/d)   + = source into aquifer,  − = sink")
    ax.set_title("CdL — mean water budget by term (post spin-up)")
    seen = {}
    for c in summ["compartment"]:
        seen.setdefault(c, COMP_COLORS.get(c, "0.6"))
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=v) for v in seen.values()],
              labels=list(seen.keys()), fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(OUT / "budget_compartment_mean.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- yearly volumes (m3/yr) ----
    days = net.index.days_in_month.to_numpy(dtype=float)
    vol = net.multiply(days, axis=0)                   # m3 per period
    yearly = vol.groupby(vol.index.year).sum()         # year-agnostic across pandas versions
    yearly.to_csv(OUT / "budget_yearly_volume_m3.csv")

    if len(yearly) >= 1:
        fig, ax = plt.subplots(figsize=(10, 5))
        bottom_pos = np.zeros(len(yearly)); bottom_neg = np.zeros(len(yearly))
        x = np.arange(len(yearly))
        for b in bases:
            vals = yearly[b].to_numpy()
            base = np.where(vals >= 0, bottom_pos, bottom_neg)
            ax.bar(x, vals, bottom=base, label=b,
                   color=COMP_COLORS.get(compartment_of(b), "0.6"),
                   edgecolor="k", lw=0.3)
            bottom_pos += np.where(vals >= 0, vals, 0)
            bottom_neg += np.where(vals < 0, vals, 0)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(yearly.index, rotation=0)
        ax.set_ylabel("volume (m³/yr)   + = source into aquifer")
        ax.set_title("CdL — yearly water budget by term")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(OUT / "budget_yearly_volume.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print(f"   wrote budget_compartment_mean.* and budget_yearly_volume.*  "
          f"({n_real} periods, {len(yearly)} yr)")


# =============================================================================
# 3b. PACKAGE-INTERNAL BUDGETS (UZF = unsaturated, SFR = surface water)
# =============================================================================
def package_budget(bud_path, title, fname):
    if not bud_path.exists():
        print(f"   !! {bud_path.name} not found. Skipping {title}.")
        return
    print(f">> [3b] {title} internal budget ({bud_path.name}) …")
    cbf = flopy.utils.CellBudgetFile(str(bud_path))
    kk = real_kstpkper(cbf) or cbf.get_kstpkper()
    raw = cbf.get_unique_record_names()
    means = {}
    for rname in raw:
        label = (rname.decode() if isinstance(rname, bytes) else rname).strip()
        if label.upper().startswith("FLOW-JA-FACE"):
            continue
        tot = []
        for k in kk:
            try:
                data = cbf.get_data(text=rname, kstpkper=k)
            except Exception:
                continue
            if not data:
                continue
            d = data[0]
            if hasattr(d, "dtype") and d.dtype.names and "q" in d.dtype.names:
                tot.append(float(np.sum(d["q"])))
            else:
                tot.append(float(np.nansum(np.asarray(d))))
        if tot:
            means[label] = np.mean(tot)
    cbf.close()
    if not means:
        print(f"   !! no records summed for {title}.")
        return
    s = pd.Series(means).sort_values()
    s.to_frame("mean_rate_m3d").to_csv(OUT / f"{fname}.csv")
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(s) + 1.5))
    ax.barh(s.index, s.values,
            color=["tab:red" if v < 0 else "tab:blue" for v in s.values],
            edgecolor="k", lw=0.4)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean rate (m³/d)   (MF6 sign: + into the package element)")
    ax.set_title(f"CdL — {title} budget (post spin-up)")
    fig.tight_layout()
    fig.savefig(OUT / f"{fname}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote {fname}.png / .csv  ({len(s)} terms)")


# =============================================================================
# 3c. PER-LAYER STORAGE CHANGE (aquifer layers) FROM THE GWF .cbc
# =============================================================================
def layer_storage():
    print(">> [3c] Per-layer storage change (GWF .cbc) …")
    if not CBC.exists():
        print(f"   !! {CBC.name} not found. Skipping.")
        return
    cbf = flopy.utils.CellBudgetFile(str(CBC))
    kk = real_kstpkper(cbf) or cbf.get_kstpkper()
    names = [(r.decode() if isinstance(r, bytes) else r).strip()
             for r in cbf.get_unique_record_names()]
    sto_terms = [n for n in names if n.upper().startswith("STO")]
    if not sto_terms:
        print("   !! no STO records (steady-state only run?). Skipping.")
        cbf.close()
        return
    per_layer = None
    for k in kk:
        tot = None
        for term in sto_terms:
            d = np.squeeze(cbf.get_data(text=term, kstpkper=k)[0])  # (nlay,ncpl)
            if d.ndim == 1:
                d = d[None, :]
            tot = d if tot is None else tot + d
        lay_sum = np.nansum(tot, axis=1)                            # (nlay,)
        per_layer = lay_sum[None, :] if per_layer is None else np.vstack([per_layer, lay_sum])
    cbf.close()
    mean_lay = np.nanmean(per_layer, axis=0)
    nlay = mean_lay.size
    pd.DataFrame({"layer": np.arange(1, nlay + 1), "mean_storage_rate_m3d": mean_lay}) \
        .to_csv(OUT / "layer_storage_mean.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(np.arange(1, nlay + 1), mean_lay, color="tab:brown", edgecolor="k")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("layer"); ax.set_ylabel("mean storage rate (m³/d)")
    ax.set_title("CdL — per-layer net storage change (+ = released to flow)")
    ax.set_xticks(np.arange(1, nlay + 1))
    fig.tight_layout()
    fig.savefig(OUT / "layer_storage_mean.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote layer_storage_mean.png / .csv  ({nlay} layers)")


# =============================================================================
# 4. LAKE (POND) STAGE + SEEPAGE
# =============================================================================
def lake_outputs():
    print(">> [4] Lake (pond) stage + seepage …")
    STG = WORKSPACE / f"{MODEL_NAME}.lak.stage"
    LKB = WORKSPACE / f"{MODEL_NAME}.lak.bud"
    if not STG.exists() and not LKB.exists():
        print("   !! no LAK output (no ponds / not run). Skipping.")
        return
    # --- stage time series per lake ---
    if STG.exists():
        try:
            sf = flopy.utils.HeadFile(str(STG), text="STAGE")
            kk = [k for k in sf.get_kstpkper() if k[1] >= SPINUP_NPER] or sf.get_kstpkper()
            stages = np.array([np.ravel(sf.get_data(kstpkper=k)) for k in kk])    # (nt, nlakes)
            stages = np.where(np.abs(stages) > 1e29, np.nan, stages)   # mask MF6 dry-lake sentinel (pond dry)
            d = real_dates(len(kk))
            nlk = stages.shape[1]
            out = pd.DataFrame(stages, columns=[f"lake{L}" for L in range(nlk)])
            out.insert(0, "date", d)
            out.to_csv(OUT / "lake_stage.csv", index=False)
            fig, ax = plt.subplots(figsize=(10, 4))
            for L in range(nlk):
                ax.plot(d, stages[:, L], lw=1.3, label=f"lake {L}")
            ax.set_ylabel("lake stage (m)"); ax.set_title("CdL — pond (LAK) stage")
            ax.legend(fontsize=8, ncol=min(nlk, 5)); ax.grid(alpha=0.3)
            fig.tight_layout(); fig.savefig(OUT / "lake_stage.png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"   wrote lake_stage.png / .csv  ({nlk} lakes)")
        except Exception as e:
            print(f"   (lake stage skipped: {e!r})")
        # --- per-pond stage panels with bottom + spill-invert datums (2026-07-04) -----------
        # metadata parsed from the WRITTEN cdl_gwf.lak (boundnames pondN, outlet inverts) and
        # the .lakN.tab files (bottom = min table stage) -> always consistent with the run.
        try:
            names, inverts, bottoms = {}, {}, {}
            blk = None
            for ln in (WORKSPACE / f"{MODEL_NAME}.lak").read_text().splitlines():
                s = ln.strip(); low = s.lower()
                if low.startswith("begin "):
                    blk = low.split()[1]; continue
                if low.startswith("end "):
                    blk = None; continue
                if not s or s.startswith("#") or blk is None:
                    continue
                p = s.split()
                if blk == "packagedata":
                    names[int(p[0]) - 1] = p[-1]
                elif blk == "outlets":
                    inverts[int(p[1]) - 1] = float(p[4])
            for L in range(nlk):
                tab = WORKSPACE / f"{MODEL_NAME}.lak{L + 1}.tab"
                if tab.exists():
                    stgs = []
                    for tl in tab.read_text().splitlines():
                        try:
                            stgs.append(float(tl.split()[0]))
                        except (ValueError, IndexError):
                            pass
                    if stgs:
                        bottoms[L] = min(stgs)
            ncols = 4
            nrows = int(np.ceil(nlk / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.1 * nrows),
                                     sharex=True, squeeze=False)
            for L in range(nlk):
                ax = axes.ravel()[L]
                ax.plot(d, stages[:, L], lw=0.7, color="tab:blue")
                if L in bottoms:
                    ax.axhline(bottoms[L], color="saddlebrown", lw=0.8, ls="--")
                if L in inverts:
                    ax.axhline(inverts[L], color="crimson", lw=0.8, ls=":")
                ttl = names.get(L, f"lake {L}")
                if L in bottoms and L in inverts:
                    ttl += f"  (bottom {bottoms[L]:.1f}, spill {inverts[L]:.1f})"
                ax.set_title(ttl, fontsize=9)
                ax.tick_params(labelsize=7)
            for ax in axes.ravel()[nlk:]:
                ax.axis("off")
            fig.suptitle("CdL — pond (LAK) stages  (brown dashes = pond bottom, red dots = spill invert)")
            fig.tight_layout()
            fig.savefig(OUT / "lake_stages_panels.png", dpi=130, bbox_inches="tight")
            plt.close(fig)
            print(f"   wrote lake_stages_panels.png ({nlk} ponds)")
        except Exception as e:
            print(f"   (lake stage panels skipped: {e!r})")
    # --- lake budget terms ('GWF' = lakebed seepage = the SW-GW flux) ---
    if LKB.exists():
        try:
            cbf = flopy.utils.CellBudgetFile(str(LKB))
            kk = [k for k in cbf.get_kstpkper() if k[1] >= SPINUP_NPER] or cbf.get_kstpkper()
            recs = [(r.decode() if isinstance(r, bytes) else r).strip()
                    for r in cbf.get_unique_record_names()]
            summary = {}
            for r in recs:
                if r.upper().startswith("FLOW-JA-FACE"):
                    continue
                tot = []
                for k in kk:
                    try:
                        data = cbf.get_data(text=r, kstpkper=k)
                    except Exception:
                        continue
                    if data and hasattr(data[0], "dtype") and data[0].dtype.names \
                            and "q" in data[0].dtype.names:
                        tot.append(float(np.sum(data[0]["q"])))
                if tot:
                    summary[r] = np.mean(tot)
            cbf.close()
            if summary:
                s = pd.Series(summary).sort_values()
                s.to_frame("mean_rate_m3d").to_csv(OUT / "lake_budget_mean.csv")
                fig, ax = plt.subplots(figsize=(8, 0.5 * len(s) + 1.5))
                ax.barh(s.index, s.values,
                        color=["tab:red" if v < 0 else "tab:blue" for v in s.values], edgecolor="k")
                ax.axvline(0, color="k", lw=0.8)
                ax.set_xlabel("mean rate (m³/d)   (MF6 sign: + into the lake)")
                ax.set_title("CdL — pond (LAK) budget  ('GWF' = bed seepage = SW–GW flux)")
                fig.tight_layout(); fig.savefig(OUT / "lake_budget_mean.png", dpi=150, bbox_inches="tight")
                plt.close(fig)
                print(f"   wrote lake_budget_mean.png / .csv  ({len(s)} terms)")
        except Exception as e:
            print(f"   (lake budget skipped: {e!r})")


# =============================================================================
if __name__ == "__main__":
    print(f">> Post-processing {MODEL_NAME} -> {OUT}")
    for fn in (obs_timeseries, head_maps, depth_maps, list_budget, layer_storage, lake_outputs):
        try:
            fn()
        except Exception as e:
            print(f"   !! {fn.__name__} failed: {e!r}")
    try:
        package_budget(UZFB, "unsaturated zone (UZF)", "uzf_budget_mean")
    except Exception as e:
        print(f"   !! uzf budget failed: {e!r}")
    try:
        package_budget(SFRB, "surface water (SFR)", "sfr_budget_mean")
    except Exception as e:
        print(f"   !! sfr budget failed: {e!r}")
    print(">> Done. Outputs in:", OUT.resolve())
