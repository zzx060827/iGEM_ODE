import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# State vector
# ---------------------------------------------------------------------
# y = [
#   A_blood,
#
#   A_liver_v,  A_liver_isf,
#   A_spleen_v, A_spleen_isf,
#   A_kidney_v, A_kidney_isf,
#   A_heart_v,  A_heart_isf,
#   A_muscle_v, A_muscle_isf,
#   A_lung_v,   A_lung_isf,
#   A_rest_v,   A_rest_isf,
#
#   B, EE, LE, CY, Ncap, Nss, Nds, Epi, M, P, Ab
# ]
#
# PBPK layer:
#   A_blood      = central blood AAV amount
#   A_{organ}_v  = vascular AAV amount in each organ
#   A_{organ}_isf = interstitial AAV amount in each organ
#
# Liver cellular layer:
#   B     = liver cell-surface bound AAV
#   EE    = early endosome
#   LE    = late endosome
#   CY    = cytosolic capsid
#   Ncap  = nuclear capsid
#   Nss   = nuclear ssDNA
#   Nds   = nuclear dsDNA
#   Epi   = episome
#   M     = mRNA
#   P     = protein
#   Ab    = simplified antibody level


ORGANS = ["liver", "spleen", "kidney", "heart", "muscle", "lung", "rest"]


def organ_fluxes(A_blood, A_v, A_isf, Ab, organ, p):
    """
    Compute PBPK fluxes for one organ.

    A_blood: central blood AAV amount
    A_v: vascular AAV amount in this organ
    A_isf: interstitial AAV amount in this organ
    Ab: antibody level
    organ: organ name string
    p: parameter dictionary

    Returns:
        J_blood_to_v
        J_v_to_isf
        J_res
        J_deg_isf
    """

    Cb = A_blood / p["V_blood"]
    Cv = A_v / p[f"V_{organ}_v"]
    Cisf = A_isf / p[f"V_{organ}_isf"]

    # Blood <-> organ vascular exchange
    J_blood_to_v = p[f"Q_{organ}"] * (Cb - Cv)

    # Organ vascular <-> organ interstitial exchange
    J_v_to_isf = p[f"PS_{organ}"] * (Cv - Cisf / p[f"Kp_{organ}"])

    # Organ vascular loss, e.g. RES uptake or nonspecific vascular clearance
    J_res = p[f"k_res_{organ}"] * A_v

    # Interstitial degradation or loss
    J_deg_isf = p[f"k_deg_isf_{organ}"] * A_isf

    return J_blood_to_v, J_v_to_isf, J_res, J_deg_isf


def rhs(t, y, p):
    # -----------------------------------------------------------------
    # Unpack state vector
    # -----------------------------------------------------------------
    (
        A_blood,

        A_liver_v,  A_liver_isf,
        A_spleen_v, A_spleen_isf,
        A_kidney_v, A_kidney_isf,
        A_heart_v,  A_heart_isf,
        A_muscle_v, A_muscle_isf,
        A_lung_v,   A_lung_isf,
        A_rest_v,   A_rest_isf,

        B, EE, LE, CY, Ncap, Nss, Nds, Epi, M, P, Ab
    ) = y

    # Store organ vascular/interstitial states in dictionaries
    A_v = {
        "liver": A_liver_v,
        "spleen": A_spleen_v,
        "kidney": A_kidney_v,
        "heart": A_heart_v,
        "muscle": A_muscle_v,
        "lung": A_lung_v,
        "rest": A_rest_v,
    }

    A_isf = {
        "liver": A_liver_isf,
        "spleen": A_spleen_isf,
        "kidney": A_kidney_isf,
        "heart": A_heart_isf,
        "muscle": A_muscle_isf,
        "lung": A_lung_isf,
        "rest": A_rest_isf,
    }

    # -----------------------------------------------------------------
    # PBPK fluxes for all organs
    # -----------------------------------------------------------------
    J_blood_to_v = {}
    J_v_to_isf = {}
    J_res = {}
    J_deg_isf = {}

    for organ in ORGANS:
        (
            J_blood_to_v[organ],
            J_v_to_isf[organ],
            J_res[organ],
            J_deg_isf[organ],
        ) = organ_fluxes(A_blood, A_v[organ], A_isf[organ], Ab, organ, p)

    # -----------------------------------------------------------------
    # Liver cell-surface binding
    # Only liver interstitial AAV drives liver cellular uptake.
    # -----------------------------------------------------------------
    C_liver_isf = A_liver_isf / p["V_liver_isf"]

    R_free = max(p["R_tot"] - B, 0.0)
    J_bind = p["k_on"] * C_liver_isf * R_free - p["k_off"] * B

    # Optional safeguard:
    # If B is numerically near zero, prevent excessive negative binding flux.
    if B <= 0 and J_bind < 0:
        J_bind = 0.0

    # -----------------------------------------------------------------
    # Antibody neutralization in central blood
    # -----------------------------------------------------------------
    J_neut_blood = p["k_neut"] * Ab * A_blood

    # -----------------------------------------------------------------
    # ODEs: central blood
    # -----------------------------------------------------------------
    # Blood loses AAV into each organ if Cb > Cv.
    # If an organ vascular concentration is higher than blood,
    # J_blood_to_v becomes negative and blood gains AAV back.
    dA_blood = (
        -sum(J_blood_to_v.values())
        - p["CL_blood"] * (A_blood / p["V_blood"])
        - J_neut_blood
    )

    # -----------------------------------------------------------------
    # ODEs: organ vascular and interstitial layers
    # -----------------------------------------------------------------
    dA_v = {}
    dA_isf = {}

    for organ in ORGANS:
        dA_v[organ] = (
            J_blood_to_v[organ]
            - J_v_to_isf[organ]
            - J_res[organ]
        )

        dA_isf[organ] = (
            J_v_to_isf[organ]
            - J_deg_isf[organ]
        )

    # Liver interstitium additionally loses AAV through productive
    # cell-surface binding.
    dA_isf["liver"] -= J_bind

    # -----------------------------------------------------------------
    # ODEs: liver intracellular fate
    # -----------------------------------------------------------------
    dB = J_bind - p["k_int"] * B

    dEE = (
        p["k_int"] * B
        - (p["k_ee_le"] + p["k_rec"] + p["k_deg_ee"]) * EE
    )

    dLE = (
        p["k_ee_le"] * EE
        - (p["k_escape"] + p["k_lys"]) * LE
    )

    dCY = (
        p["k_escape"] * LE
        - (p["k_nuc"] + p["k_uncoat_cyto"] + p["k_deg_cyto"]) * CY
    )

    dNcap = (
        p["k_nuc"] * CY
        - (p["k_uncoat_nuc"] + p["k_deg_ncap"]) * Ncap
    )

    dNss = (
        p["k_uncoat_cyto"] * CY
        + p["k_uncoat_nuc"] * Ncap
        - (p["k_ds"] + p["k_deg_ss"]) * Nss
    )

    dNds = (
        p["k_ds"] * Nss
        - (p["k_epi"] + p["k_deg_ds"]) * Nds
    )

    dEpi = (
        p["k_epi"] * Nds
        - (p["k_loss_epi"] + p["k_dil"]) * Epi
    )

    # -----------------------------------------------------------------
    # Expression module
    # -----------------------------------------------------------------
    tx = (
        p["k_tx"]
        * (Epi ** p["h"])
        / (p["EC50_tx"] ** p["h"] + Epi ** p["h"] + 1e-30)
    )

    dM = tx - p["k_deg_m"] * M
    dP = p["k_tl"] * M - p["k_deg_p"] * P

    # -----------------------------------------------------------------
    # Simplified antibody induction
    # -----------------------------------------------------------------
    # Antigen signal can include blood AAV, organ vascular AAV,
    # organ interstitial AAV, and intracellular late endosome AAV.
    Ag = (
        A_blood
        + sum(A_v.values())
        + 0.5 * sum(A_isf.values())
        + LE
    )

    dAb = (
        p["k_ab_max"] * Ag / (p["EC50_ab"] + Ag + 1e-30)
        - p["k_deg_ab"] * Ab
    )

    # -----------------------------------------------------------------
    # Return derivatives in the same order as y
    # -----------------------------------------------------------------
    return [
        dA_blood,

        dA_v["liver"],  dA_isf["liver"],
        dA_v["spleen"], dA_isf["spleen"],
        dA_v["kidney"], dA_isf["kidney"],
        dA_v["heart"],  dA_isf["heart"],
        dA_v["muscle"], dA_isf["muscle"],
        dA_v["lung"],   dA_isf["lung"],
        dA_v["rest"],   dA_isf["rest"],

        dB, dEE, dLE, dCY, dNcap, dNss, dNds, dEpi, dM, dP, dAb
    ]


# ---------------------------------------------------------------------
# Parameters(simulated in a 25 g mouse)
# ---------------------------------------------------------------------
p = {
    # -------------------------------------------------------------
    # Central blood
    # -------------------------------------------------------------
    "V_blood": 1.5,      # mL, total amont of blood in the mouse
    "CL_blood": 0.01,    # mL/h-like effective clearance

    # -------------------------------------------------------------
    # Organ volumes
    # Vascular and interstitial effective volumes
    # -------------------------------------------------------------

# considering the vascular volume and interstitial fluid volume take 10% and 20% of the total organ volume.
    "V_liver_v": 0.14,
    "V_liver_isf": 0.22,

    "V_spleen_v": 0.01,
    "V_spleen_isf": 0.02,

    "V_kidney_v": 0.05,
    "V_kidney_isf": 0.08,

    "V_heart_v": 0.012,
    "V_heart_isf": 0.025,

    "V_muscle_v": 0.4,
    "V_muscle_isf": 1.2,

    "V_lung_v": 0.08,
    "V_lung_isf": 0.04,

    "V_rest_v": 0.7,
    "V_rest_isf": 2.16,

    # -------------------------------------------------------------
    # Organ blood-flow-like exchange terms
    # Units should be compatible with mL/h.
    # -------------------------------------------------------------
# considering the total blood flow is 500 ml/h
    'CO': 500,  # mL/h, cardiac output for a 25 g mouse
    "Q_lung": 500.0,  

    "Q_liver": 0.25 * 500.0,     # 肝脏约占 CO 的 25%
    "Q_spleen": 0.06 * 500.0,    # 脾脏约占 CO 的 6%
    "Q_kidney": 0.20 * 500.0,    # 肾脏约占 CO 的 20%
    "Q_heart": 0.05 * 500.0,     # 心脏自身血流约 5%
    "Q_muscle": 0.15 * 500.0,    # 骨骼肌约占 CO 的 15%
    "Q_rest": 0.14 * 500.0,       # 剩余小器官约 14%

    # -------------------------------------------------------------
    # Vascular-to-interstitial permeability terms
    # Higher PS means easier extravasation into tissue interstitium.
    # -------------------------------------------------------------
    "PS_liver": 0.20,
    "PS_spleen": 0.15,
    "PS_kidney": 0.08,
    "PS_heart": 0.05,
    "PS_muscle": 0.03,
    "PS_lung": 0.10,
    "PS_rest": 0.02,

    # -------------------------------------------------------------
    # Tissue partition coefficients
    # Higher Kp means more retention in interstitial/tissue side.
    # -------------------------------------------------------------
    "Kp_liver": 1.5,
    "Kp_spleen": 1.2,
    "Kp_kidney": 0.8,
    "Kp_heart": 0.6,
    "Kp_muscle": 0.5,
    "Kp_lung": 0.7,
    "Kp_rest": 0.4,

    # -------------------------------------------------------------
    # Organ vascular RES/nonspecific clearance
    # Liver and spleen are set higher because they often dominate
    # reticuloendothelial uptake.
    # -------------------------------------------------------------
    "k_res_liver": 0.02,
    "k_res_spleen": 0.03,
    "k_res_kidney": 0.005,
    "k_res_heart": 0.002,
    "k_res_muscle": 0.001,
    "k_res_lung": 0.004,
    "k_res_rest": 0.001,

    # -------------------------------------------------------------
    # Interstitial degradation/loss
    # -------------------------------------------------------------
    "k_deg_isf_liver": 0.001,
    "k_deg_isf_spleen": 0.001,
    "k_deg_isf_kidney": 0.001,
    "k_deg_isf_heart": 0.001,
    "k_deg_isf_muscle": 0.001,
    "k_deg_isf_lung": 0.001,
    "k_deg_isf_rest": 0.001,

    # ------------------------------------------------------------- 
    # Liver cell-surface binding
    # -------------------------------------------------------------
    "R_tot": 1e5,        # a.u. or copies/cell-equivalent
    "k_on": 1e-6,
    "k_off": 0.05,
    "k_int": 0.2,

    # -------------------------------------------------------------
    # Liver intracellular trafficking
    # -------------------------------------------------------------
    "k_ee_le": 0.3,
    "k_rec": 0.05,
    "k_deg_ee": 0.02,
    "k_escape": 0.005,
    "k_lys": 0.1,
    "k_nuc": 0.02,
    "k_uncoat_cyto": 0.005,
    "k_uncoat_nuc": 0.02,
    "k_deg_cyto": 0.01,
    "k_deg_ncap": 0.005,
    "k_ds": 0.01,
    "k_deg_ss": 0.02,
    "k_epi": 0.01,
    "k_deg_ds": 0.005,
    "k_loss_epi": 0.0005,
    "k_dil": 0.0,

    # -------------------------------------------------------------
    # Expression
    # -------------------------------------------------------------
    "k_tx": 2.0,
    "h": 1.2,
    "EC50_tx": 100.0,
    "k_deg_m": 0.2,
    "k_tl": 5.0,
    "k_deg_p": 0.02,

    # -------------------------------------------------------------
    # Antibody
    # -------------------------------------------------------------
    "k_neut": 1e-14,
    "k_ab_max": 0.05,
    "EC50_ab": 1e10,
    "k_deg_ab": 0.005,
}


# ---------------------------------------------------------------------
# Initial condition
# ---------------------------------------------------------------------
dose_vg = 1e12

y0 = [
    dose_vg,

    0, 0,   # liver_v, liver_isf
    0, 0,   # spleen_v, spleen_isf
    0, 0,   # kidney_v, kidney_isf
    0, 0,   # heart_v, heart_isf
    0, 0,   # muscle_v, muscle_isf
    0, 0,   # lung_v, lung_isf
    0, 0,   # rest_v, rest_isf

    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    # B, EE, LE, CY, Ncap, Nss, Nds, Epi, M, P, Ab
]


# ---------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------

IDX = {
    "A_blood": 0,

    "A_liver_v": 1,
    "A_liver_isf": 2,
    "A_spleen_v": 3,
    "A_spleen_isf": 4,
    "A_kidney_v": 5,
    "A_kidney_isf": 6,
    "A_heart_v": 7,
    "A_heart_isf": 8,
    "A_muscle_v": 9,
    "A_muscle_isf": 10,
    "A_lung_v": 11,
    "A_lung_isf": 12,
    "A_rest_v": 13,
    "A_rest_isf": 14,

    "B": 15,
    "EE": 16,
    "LE": 17,
    "CY": 18,
    "Ncap": 19,
    "Nss": 20,
    "Nds": 21,
    "Epi": 22,
    "M": 23,
    "P": 24,
    "Ab": 25,
}


## ---------------------------------------------------------------------
# Solve: 短时间全身分布 + 长时间胞内/表达模拟
# ---------------------------------------------------------------------

# 短时间：用于看血液、器官血管、间质快速分布
t_eval_short = np.linspace(0, 2, 300)  # 0-2 h

sol_short = solve_ivp(
    lambda t, y: rhs(t, y, p),
    t_span=(0, t_eval_short[-1]),
    y0=y0,
    t_eval=t_eval_short,
    method="Radau",
    rtol=1e-6,
    atol=1e-9,
)

if not sol_short.success:
    raise RuntimeError(sol_short.message)


# 长时间：用于看肝细胞内 trafficking、episome、mRNA、protein、antibody
t_eval_long = np.linspace(0, 24 * 56, 800)  # 56 days, hours

sol_long = solve_ivp(
    lambda t, y: rhs(t, y, p),
    t_span=(0, t_eval_long[-1]),
    y0=y0,
    t_eval=t_eval_long,
    method="Radau",
    rtol=1e-6,
    atol=1e-9,
)

if not sol_long.success:
    raise RuntimeError(sol_long.message)


# ---------------------------------------------------------------------
# 绘图设置
# ---------------------------------------------------------------------

colors = {
    "liver": "red",
    "spleen": "green",
    "kidney": "blue",
    "heart": "orange",
    "muscle": "purple",
    "lung": "cyan",
    "rest": "brown"
}

def positive_only(x, floor=1.0):
    """
    log-scale 绘图用。
    把 <=0 的值替换成 floor，避免 log(0) 报错。
    """
    x = np.asarray(x)
    return np.where(x > floor, x, floor)


# ---------------------------------------------------------------------
# 1. 短时间全身血管/间质浓度：0-2 h
# ---------------------------------------------------------------------

plt.figure(figsize=(14, 6))

# -----------------------------
# 血管浓度
# -----------------------------
plt.subplot(1, 2, 1)

for organ in ORGANS:
    idx_v = IDX[f"A_{organ}_v"]
    C_v = sol_short.y[idx_v] / p[f"V_{organ}_v"]

    plt.plot(
        sol_short.t,
        positive_only(C_v),
        label=f"{organ}_vascular",
        color=colors[organ],
        linestyle="-",
        linewidth=1.8
    )

C_blood_short = sol_short.y[IDX["A_blood"]] / p["V_blood"]

plt.plot(
    sol_short.t,
    positive_only(C_blood_short),
    label="blood",
    color="black",
    linestyle="--",
    linewidth=2.0
)

plt.yscale("log")
plt.xlabel("Time (h)")
plt.ylabel("AAV concentration (vg/mL)")
plt.title("Vascular AAV concentrations, 0-2 h")
plt.grid(True, which="both", linestyle="--", alpha=0.35)
plt.legend(fontsize=8)


# -----------------------------
# 间质浓度
# -----------------------------
plt.subplot(1, 2, 2)

for organ in ORGANS:
    idx_isf = IDX[f"A_{organ}_isf"]
    C_isf = sol_short.y[idx_isf] / p[f"V_{organ}_isf"]

    plt.plot(
        sol_short.t,
        positive_only(C_isf),
        label=f"{organ}_ISF",
        color=colors[organ],
        linewidth=1.8
    )

plt.yscale("log")
plt.xlabel("Time (h)")
plt.ylabel("AAV concentration (vg/mL)")
plt.title("Interstitial AAV concentrations, 0-2 h")
plt.grid(True, which="both", linestyle="--", alpha=0.35)
plt.legend(fontsize=8)

plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------
# 2. 长时间肝细胞内 trafficking 各层变化：0-56 days
# ---------------------------------------------------------------------

t_day_long = sol_long.t / 24.0

plt.figure(figsize=(13, 6))

intracellular_states = [
    ("B", "Surface-bound AAV"),
    ("EE", "Early endosome"),
    ("LE", "Late endosome"),
    ("CY", "Cytosolic capsid"),
    ("Ncap", "Nuclear capsid"),
    ("Nss", "Nuclear ssDNA"),
    ("Nds", "Nuclear dsDNA"),
    ("Epi", "Episome"),
]

for state, label in intracellular_states:
    plt.plot(
        t_day_long,
        positive_only(sol_long.y[IDX[state]]),
        label=label,
        linewidth=1.8
    )

plt.yscale("log")
plt.xlabel("Time (day)")
plt.ylabel("Amount / arbitrary units")
plt.title("Liver intracellular AAV trafficking states, 0-56 days")
plt.grid(True, which="both", linestyle="--", alpha=0.35)
plt.legend(fontsize=8)
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------
# 3. 长时间表达模块：Epi, mRNA, Protein
# ---------------------------------------------------------------------

plt.figure(figsize=(12, 6))

plt.plot(
    t_day_long,
    positive_only(sol_long.y[IDX["Epi"]]),
    label="Episome (Epi)",
    color="red",
    linewidth=2.0
)

plt.plot(
    t_day_long,
    positive_only(sol_long.y[IDX["M"]]),
    label="mRNA (M)",
    color="blue",
    linewidth=2.0
)

plt.plot(
    t_day_long,
    positive_only(sol_long.y[IDX["P"]]),
    label="Protein (P)",
    color="green",
    linewidth=2.0
)

plt.yscale("log")
plt.xlabel("Time (day)")
plt.ylabel("Amount / arbitrary units")
plt.title("Liver expression states, 0-56 days")
plt.grid(True, which="both", linestyle="--", alpha=0.35)
plt.legend()
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------
# 4. 长时间抗体水平：0-56 days
# ---------------------------------------------------------------------

plt.figure(figsize=(8, 4))

plt.plot(
    t_day_long,
    positive_only(sol_long.y[IDX["Ab"]]),
    color="magenta",
    linewidth=2.0
)

plt.yscale("log")
plt.xlabel("Time (day)")
plt.ylabel("Antibody level (Ab, a.u.)")
plt.title("Simplified antibody kinetics, 0-56 days")
plt.grid(True, which="both", linestyle="--", alpha=0.35)
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------
# 5. 输出指标：用长时间结果计算
# ---------------------------------------------------------------------

C_blood_long = sol_long.y[IDX["A_blood"]] / p["V_blood"]

AUC_blood = np.trapezoid(C_blood_long, sol_long.t)
Cmax_blood = np.max(C_blood_long)

peak_liver_isf = np.max(sol_long.y[IDX["A_liver_isf"]])
peak_spleen_isf = np.max(sol_long.y[IDX["A_spleen_isf"]])
peak_kidney_isf = np.max(sol_long.y[IDX["A_kidney_isf"]])
peak_heart_isf = np.max(sol_long.y[IDX["A_heart_isf"]])
peak_muscle_isf = np.max(sol_long.y[IDX["A_muscle_isf"]])
peak_lung_isf = np.max(sol_long.y[IDX["A_lung_isf"]])
peak_rest_isf = np.max(sol_long.y[IDX["A_rest_isf"]])

peak_epi = np.max(sol_long.y[IDX["Epi"]])
peak_mRNA = np.max(sol_long.y[IDX["M"]])
peak_protein = np.max(sol_long.y[IDX["P"]])
peak_antibody = np.max(sol_long.y[IDX["Ab"]])

print("AUC_blood:", AUC_blood)
print("Cmax_blood:", Cmax_blood)
print("Peak liver ISF:", peak_liver_isf)
print("Peak spleen ISF:", peak_spleen_isf)
print("Peak kidney ISF:", peak_kidney_isf)
print("Peak heart ISF:", peak_heart_isf)
print("Peak muscle ISF:", peak_muscle_isf)
print("Peak lung ISF:", peak_lung_isf)
print("Peak rest ISF:", peak_rest_isf)
print("Peak liver episome:", peak_epi)
print("Peak mRNA:", peak_mRNA)
print("Peak protein:", peak_protein)
print("Peak antibody:", peak_antibody)