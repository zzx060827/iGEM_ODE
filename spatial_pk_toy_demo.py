"""
Minimal 1D spatial-PK demo for AAV delivery.

This file is intentionally dependency-free: it uses only the Python standard
library and writes CSV + SVG outputs. It is a bridge between the current 0D PBPK
ODE model and a future CFD / advection-diffusion-reaction spatial model.

Model:
    vascular concentration C(x,t)
    receptor-bound AAV B(x,t)
    internalized AAV I(x,t)
    episomal / expression-competent vector E(x,t)

    dC/dt = -u dC/dx + D d2C/dx2 - uptake - nonspecific loss
    dB/dt = kon C (Bmax - B) - koff B - kint B
    dI/dt = kint B - (kescape + klys) I
    dE/dt = kescape I - kloss E

Three scenarios are compared:
    baseline_fast_iv: short systemic pulse, faster flow, lower wall access.
    capsid_enhanced_same_input: same input/flow, higher receptor wall access.
    slow_flow_local_trapping: longer local input and slower flow; useful for
        showing that residence-time control can create proximal trapping.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


OUTDIR = Path("spatial_pk_demo_outputs")


def inlet_profile(t_h: float, dose_scale: float, tinf_h: float, half_life_h: float) -> float:
    if 0.0 <= t_h <= tinf_h:
        return dose_scale / tinf_h
    return (dose_scale / tinf_h) * math.exp(-math.log(2.0) * (t_h - tinf_h) / half_life_h)


def simulate(label: str, params: dict[str, float]) -> dict[str, object]:
    n = int(params["n"])
    length_cm = params["length_cm"]
    dx = length_cm / (n - 1)
    t_end_h = params["t_end_h"]
    dt_h = params["dt_h"]
    steps = int(t_end_h / dt_h)

    c = [0.0] * n
    b = [0.0] * n
    internal = [0.0] * n
    epi = [0.0] * n

    snapshots: list[tuple[float, list[float], list[float]]] = []
    snapshot_every = max(1, steps // 120)

    for step in range(steps + 1):
        t_h = step * dt_h
        c_in = inlet_profile(t_h, params["dose_scale"], params["tinf_h"], params["input_half_life_h"])

        if step % snapshot_every == 0:
            snapshots.append((t_h, c[:], epi[:]))

        old_c = c[:]
        old_b = b[:]
        old_i = internal[:]
        old_e = epi[:]

        c[0] = c_in
        for j in range(1, n - 1):
            adv = -params["u_cm_h"] * (old_c[j] - old_c[j - 1]) / dx
            diff = params["D_cm2_h"] * (old_c[j + 1] - 2.0 * old_c[j] + old_c[j - 1]) / (dx * dx)
            receptor_free = max(params["Bmax"] - old_b[j], 0.0)
            bind = params["kon"] * old_c[j] * receptor_free
            unbind = params["koff"] * old_b[j]
            uptake = params["wall_access"] * max(bind - unbind, 0.0)
            loss = params["k_loss_extra"] * old_c[j]
            c[j] = max(old_c[j] + dt_h * (adv + diff - uptake - loss), 0.0)

        c[-1] = c[-2]

        for j in range(n):
            receptor_free = max(params["Bmax"] - old_b[j], 0.0)
            bind = params["wall_access"] * params["kon"] * old_c[j] * receptor_free
            unbind = params["koff"] * old_b[j]
            internalize = params["kint"] * old_b[j]
            escape = params["kescape"] * old_i[j]
            lys = params["klys"] * old_i[j]

            b[j] = max(old_b[j] + dt_h * (bind - unbind - internalize), 0.0)
            internal[j] = max(old_i[j] + dt_h * (internalize - escape - lys), 0.0)
            epi[j] = max(old_e[j] + dt_h * (escape - params["kloss_epi"] * old_e[j]), 0.0)

    mean_epi = sum(epi) / n
    sd_epi = math.sqrt(sum((x - mean_epi) ** 2 for x in epi) / n)
    distal = sum(epi[int(0.75 * n):]) / max(sum(epi), 1e-30)
    return {
        "label": label,
        "x": [j * dx for j in range(n)],
        "c": c,
        "b": b,
        "internal": internal,
        "epi": epi,
        "snapshots": snapshots,
        "metrics": {
            "total_epi": sum(epi),
            "mean_epi": mean_epi,
            "cv_epi": sd_epi / max(mean_epi, 1e-30),
            "distal_quarter_epi_fraction": distal,
            "peak_epi": max(epi),
        },
    }


def write_profile_csv(result: dict[str, object]) -> None:
    OUTDIR.mkdir(exist_ok=True)
    path = OUTDIR / f"{result['label']}_final_profile.csv"
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x_cm", "C_final", "B_bound", "I_internal", "E_epi"])
        for row in zip(result["x"], result["c"], result["b"], result["internal"], result["epi"]):
            writer.writerow([f"{v:.8g}" for v in row])


def color(value: float, vmin: float, vmax: float) -> str:
    z = 0.0 if vmax <= vmin else (value - vmin) / (vmax - vmin)
    z = max(0.0, min(1.0, z))
    r = int(255 * z)
    g = int(70 + 90 * (1.0 - abs(z - 0.5) * 2.0))
    b = int(255 * (1.0 - z))
    return f"rgb({r},{g},{b})"


def write_svg(results: list[dict[str, object]]) -> None:
    OUTDIR.mkdir(exist_ok=True)
    width = 1100
    height = 950
    margin = 70
    panel_h = 230
    gap = 55

    all_epi = [v for r in results for v in r["epi"]]
    vmax = max(all_epi) if all_epi else 1.0

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="70" y="34" font-family="Arial" font-size="24" font-weight="700">1D AAV Spatial PK demo: flow changes spatial expression</text>',
        '<text x="70" y="58" font-family="Arial" font-size="13" fill="#444">Color shows final episomal/expression-competent vector E(x). Lower CV means more spatially uniform expression.</text>',
    ]

    for idx, result in enumerate(results):
        top = margin + idx * (panel_h + gap)
        xvals = result["x"]
        epi = result["epi"]
        metrics = result["metrics"]
        label = result["label"]
        plot_w = width - 2 * margin
        plot_h = 120
        y0 = top + 70
        dx_px = plot_w / len(xvals)

        parts.append(f'<text x="{margin}" y="{top + 20}" font-family="Arial" font-size="18" font-weight="700">{label}</text>')
        parts.append(
            f'<text x="{margin}" y="{top + 43}" font-family="Arial" font-size="13" fill="#333">'
            f'total E={metrics["total_epi"]:.3g}, CV={metrics["cv_epi"]:.3g}, distal-quarter fraction={metrics["distal_quarter_epi_fraction"]:.3g}'
            '</text>'
        )
        parts.append(f'<rect x="{margin}" y="{y0}" width="{plot_w}" height="{plot_h}" fill="#f4f4f4" stroke="#222" stroke-width="1"/>')
        for j, val in enumerate(epi):
            x = margin + j * dx_px
            parts.append(f'<rect x="{x:.2f}" y="{y0}" width="{dx_px + 0.5:.2f}" height="{plot_h}" fill="{color(val, 0.0, vmax)}"/>')

        # Line trace over heat strip.
        max_epi = max(epi) if epi else 1.0
        points = []
        for j, val in enumerate(epi):
            x = margin + j * dx_px
            y = y0 + plot_h - (val / max(max_epi, 1e-30)) * (plot_h - 8) - 4
            points.append(f"{x:.2f},{y:.2f}")
        parts.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="#111" stroke-width="2"/>')
        parts.append(f'<text x="{margin}" y="{y0 + plot_h + 28}" font-family="Arial" font-size="12">inlet / proximal</text>')
        parts.append(f'<text x="{width - margin - 90}" y="{y0 + plot_h + 28}" font-family="Arial" font-size="12">outlet / distal</text>')

    # Legend
    lx, ly = margin, height - 45
    parts.append(f'<text x="{lx}" y="{ly - 8}" font-family="Arial" font-size="12" fill="#333">low E</text>')
    for k in range(80):
        parts.append(f'<rect x="{lx + 48 + k * 3}" y="{ly - 20}" width="3" height="16" fill="{color(k, 0, 79)}"/>')
    parts.append(f'<text x="{lx + 300}" y="{ly - 8}" font-family="Arial" font-size="12" fill="#333">high E</text>')
    parts.append("</svg>")
    (OUTDIR / "spatial_pk_comparison.svg").write_text("\n".join(parts))


def write_metrics(results: list[dict[str, object]]) -> None:
    OUTDIR.mkdir(exist_ok=True)
    with (OUTDIR / "scenario_metrics.csv").open("w", newline="") as f:
        fields = ["scenario", "total_epi", "mean_epi", "cv_epi", "distal_quarter_epi_fraction", "peak_epi"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for result in results:
            row = {"scenario": result["label"]}
            row.update(result["metrics"])
            writer.writerow(row)


def main() -> None:
    base = {
        "n": 180,
        "length_cm": 1.0,
        "t_end_h": 6.0,
        "dt_h": 0.0007,
        "D_cm2_h": 0.00045,
        "dose_scale": 1.0,
        "input_half_life_h": 1.6,
        "Bmax": 1.0,
        "kon": 0.9,
        "koff": 0.12,
        "kint": 0.45,
        "kescape": 0.045,
        "klys": 0.18,
        "kloss_epi": 0.018,
        "k_loss_extra": 0.015,
    }
    scenarios = [
        ("baseline_fast_iv", {**base, "u_cm_h": 2.6, "tinf_h": 0.16, "wall_access": 0.55}),
        ("capsid_enhanced_same_input", {**base, "u_cm_h": 2.6, "tinf_h": 0.16, "wall_access": 1.20}),
        ("slow_flow_local_trapping", {**base, "u_cm_h": 1.1, "tinf_h": 0.65, "wall_access": 0.90}),
    ]
    results = [simulate(label, params) for label, params in scenarios]
    for result in results:
        write_profile_csv(result)
    write_metrics(results)
    write_svg(results)

    print("Spatial PK demo written to:", OUTDIR.resolve())
    for result in results:
        print(result["label"], result["metrics"])


if __name__ == "__main__":
    main()
