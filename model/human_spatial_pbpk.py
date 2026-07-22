"""Equation-aligned reference-human multiregion AAV PBPK model.

The extracellular and intracellular equations follow ``ode1.0.py``: effective
flow-limited vascular exchange, PS/Kp-limited ISF exchange, organ-specific
vascular/ISF loss, receptor uptake, endosomal trafficking, nuclear processing,
episome formation, and mRNA/protein turnover. Human physiology changes the
volumes, flows, dose, organ subdivisions, and CNS depth parameters, not the
underlying calculation method.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.integrate import solve_ivp


BODY_WEIGHT_KG = 70.0
DOSE_VG_PER_KG = 4.0e13
DOSE_VG = BODY_WEIGHT_KG * DOSE_VG_PER_KG
INFUSION_DURATION_H = 10.0 / 60.0
CARDIAC_OUTPUT_ML_H = 330_000.0
# Match the dimensionless Q scaling used by the mouse model. Human cardiac
# output and compartment volumes still set the species-specific transit times.
EFFECTIVE_FLOW_SCALE = 0.05
BLOOD_CAPSID_HALF_LIFE_H = 2.4

CAPSID_HALF_LIFE_H = {
    "brain": 24.8,
    "heart": 15.6,
    "liver": 22.6,
    "spleen": 22.9,
    "kidney": 24.0,
    "muscle": 48.7,
    "lung": 18.0,
    "rest": 34.0,
}

RECEPTOR_DENSITY_VG_ML = {
    "brain": 1.0e7,
    "heart": 7.0e7,
    "liver": 2.0e8,
    "spleen": 1.6e8,
    "kidney": 1.0e8,
    "muscle": 4.0e7,
    "lung": 8.0e7,
    "rest": 2.0e7,
}

BINDING_ON_ML_VG_H = {
    "brain": 5.0e-12,
    "heart": 7.0e-12,
    "liver": 1.2e-11,
    "spleen": 1.0e-11,
    "kidney": 8.0e-12,
    "muscle": 5.0e-12,
    "lung": 8.0e-12,
    "rest": 4.0e-12,
}


@dataclass(frozen=True)
class Region:
    label: str
    parent: str
    flow_fraction: float
    vascular_ml: float
    isf_ml: float
    ps_ml_h: float
    kp: float
    internalization_half_life_h: float
    episome_half_life_days: float
    route: str = "systemic"


@dataclass(frozen=True)
class AdministrationRoute:
    label: str
    label_zh: str
    infusion_duration_h: float
    input_state: str
    route_class: str
    evidence_source: str
    description: str
    description_zh: str


ADMINISTRATION_ROUTES: dict[str, AdministrationRoute] = {
    "iv": AdministrationRoute(
        label="Peripheral intravenous",
        label_zh="外周静脉注射",
        infusion_duration_h=10.0 / 60.0,
        input_state="A_arm_vein",
        route_class="systemic",
        evidence_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/",
        description="Left arm vein to right heart, lung, left heart, then systemic organs",
        description_zh="左臂静脉 → 右心 → 肺 → 左心 → 全身器官",
    ),
    "intrathecal": AdministrationRoute(
        label="Lumbar intrathecal",
        label_zh="腰椎鞘内注射",
        infusion_duration_h=10.0 / 60.0,
        input_state="A_csf_lumbar",
        route_class="csf",
        evidence_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC3618620/",
        description="Lumbar CSF to spinal/cranial CSF and CNS ISF, with venous CSF drainage",
        description_zh="腰段 CSF → 脊髓/颅内 CSF → CNS ISF，并经 CSF 回流进入静脉",
    ),
    "intramuscular": AdministrationRoute(
        label="Deltoid intramuscular",
        label_zh="三角肌肌内注射",
        infusion_duration_h=5.0 / 60.0,
        input_state="A_im_depot",
        route_class="local_depot",
        evidence_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC4098646/",
        description="Injected deltoid depot to local muscle ISF with slower systemic escape",
        description_zh="注射侧三角肌 depot → 局部肌肉 ISF，并缓慢进入全身循环",
    ),
    "intracisternal": AdministrationRoute(
        label="Intra-cisterna magna",
        label_zh="枕大池注射（ICM）",
        infusion_duration_h=10.0 / 60.0,
        input_state="A_csf_cranial",
        route_class="csf",
        evidence_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC12509745/",
        description="Cisterna magna to cranial CSF, favoring cerebellar and brainstem surfaces",
        description_zh="枕大池 → 颅内 CSF，优先接触小脑与脑干表面",
    ),
    "intracerebroventricular": AdministrationRoute(
        label="Intracerebroventricular",
        label_zh="脑室内注射（ICV）",
        infusion_duration_h=10.0 / 60.0,
        input_state="A_csf_cranial",
        route_class="csf",
        evidence_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC8417503/",
        description="Ventricular CSF input with greater access to periventricular and deep-gray regions",
        description_zh="脑室 CSF 起始，提高脑室周围与深部灰质入口",
    ),
    "inhaled": AdministrationRoute(
        label="Inhaled / airway",
        label_zh="吸入 / 气道给药",
        infusion_duration_h=15.0 / 60.0,
        input_state="A_airway_depot",
        route_class="local_depot",
        evidence_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC5655841/",
        description="Airway depot to lung ISF with limited leakage into pulmonary venous blood",
        description_zh="气道 depot → 肺部 ISF，并有少量进入肺静脉血",
    ),
}

# Adult reference CSF volumes and deliberately coarse transfer priors. These
# are mechanistic placeholders for fitting, not clinical route calibration.
ROUTE_COMPARTMENT_VOLUMES_ML = {
    "csf_lumbar": 75.0,
    "csf_cranial": 75.0,
}
CSF_ROSTRAL_FLOW_ML_H = 10.0
CSF_ABSORPTION_HALF_LIFE_H = 5.5
IM_DEPOT_ABSORPTION_HALF_LIFE_H = 3.0
IM_LOCAL_FRACTION = 0.85
AIRWAY_DEPOT_RELEASE_HALF_LIFE_H = 1.2
AIRWAY_LOCAL_FRACTION = 0.90

# CSF-to-ISF permeability and partition priors distinguish surface-accessible
# cortex/spinal cord from deep parenchyma.
CSF_REGION_ACCESS = {
    "brain_frontal": ("cranial", 0.14, 0.12),
    "brain_parietal": ("cranial", 0.12, 0.11),
    "brain_temporal": ("cranial", 0.09, 0.09),
    "brain_occipital": ("cranial", 0.08, 0.08),
    "brain_deep_gray": ("cranial", 0.018, 0.025),
    "brain_cerebellum": ("cranial", 0.10, 0.10),
    "brainstem_spinal": ("lumbar", 0.24, 0.16),
}

# Route-specific multipliers preserve one CSF transport equation while changing
# the injection origin. They are mechanistic priors that must be fitted to
# route-resolved large-animal biodistribution data.
CSF_ROUTE_ACCESS_MULTIPLIER = {
    "intrathecal": {},
    "intracisternal": {
        "brain_frontal": 1.15,
        "brain_parietal": 1.20,
        "brain_temporal": 1.10,
        "brain_occipital": 1.25,
        "brain_deep_gray": 0.70,
        "brain_cerebellum": 1.90,
        "brainstem_spinal": 1.60,
    },
    "intracerebroventricular": {
        "brain_frontal": 1.25,
        "brain_parietal": 1.10,
        "brain_temporal": 1.05,
        "brain_occipital": 0.90,
        "brain_deep_gray": 2.20,
        "brain_cerebellum": 0.65,
        "brainstem_spinal": 0.50,
    },
}


def _region(
    label: str,
    parent: str,
    flow_fraction: float,
    vascular_ml: float,
    isf_ml: float,
    ps_density: float,
    kp: float,
    internalization_half_life_h: float,
    episome_half_life_days: float,
    route: str = "systemic",
) -> Region:
    return Region(
        label=label,
        parent=parent,
        flow_fraction=flow_fraction,
        vascular_ml=vascular_ml,
        isf_ml=isf_ml,
        ps_ml_h=ps_density * isf_ml / 1000.0,
        kp=kp,
        internalization_half_life_h=internalization_half_life_h,
        episome_half_life_days=episome_half_life_days,
        route=route,
    )


# Systemic flow fractions sum to one cardiac output. Gut and spleen drain via
# the portal vein into liver; ``liver`` therefore receives 0.06 CO directly and
# 0.17 CO through the portal circulation.
REGIONS: dict[str, Region] = {
    "brain_frontal": _region("Frontal cortex", "brain", 0.025, 35, 210, 0.40, 0.050, 42, 365),
    "brain_parietal": _region("Parietal cortex", "brain", 0.024, 34, 195, 0.36, 0.048, 44, 365),
    "brain_temporal": _region("Temporal cortex", "brain", 0.021, 31, 180, 0.32, 0.045, 46, 365),
    "brain_occipital": _region("Occipital cortex", "brain", 0.018, 27, 145, 0.30, 0.043, 48, 365),
    "brain_deep_gray": _region("Deep gray nuclei", "brain", 0.018, 27, 150, 0.12, 0.028, 60, 365),
    "brain_cerebellum": _region("Cerebellum", "brain", 0.014, 22, 125, 0.18, 0.035, 54, 365),
    "brainstem_spinal": _region("Brainstem and spinal cord", "brain", 0.010, 18, 120, 0.10, 0.024, 72, 365),
    "heart": _region("Heart", "heart", 0.040, 110, 230, 20.0, 0.70, 24, 300),
    "liver": _region("Liver", "liver", 0.060, 420, 1050, 45.0, 1.60, 8, 120),
    "spleen": _region("Spleen", "spleen", 0.030, 45, 115, 35.0, 1.25, 8, 150),
    "kidney_left_cortex": _region("Left kidney cortex", "kidney", 0.060, 55, 85, 30.0, 0.85, 18, 90),
    "kidney_left_medulla": _region("Left kidney medulla", "kidney", 0.035, 32, 58, 18.0, 0.70, 24, 90),
    "kidney_right_cortex": _region("Right kidney cortex", "kidney", 0.060, 55, 85, 30.0, 0.85, 18, 90),
    "kidney_right_medulla": _region("Right kidney medulla", "kidney", 0.035, 32, 58, 18.0, 0.70, 24, 90),
    "muscle_injected_arm": _region("Injected-side arm muscle", "muscle", 0.012, 38, 420, 8.0, 0.55, 30, 365),
    "muscle_contralateral_arm": _region("Contralateral arm muscle", "muscle", 0.012, 38, 420, 8.0, 0.55, 36, 365),
    "muscle_trunk": _region("Trunk muscle", "muscle", 0.056, 175, 2900, 7.0, 0.55, 36, 365),
    "muscle_legs": _region("Leg muscle", "muscle", 0.090, 280, 5000, 7.0, 0.55, 36, 365),
    "gut": _region("Gastrointestinal tract", "rest", 0.140, 310, 1050, 18.0, 0.65, 20, 120),
    "skin_adipose": _region("Skin and adipose", "rest", 0.080, 230, 7200, 4.0, 0.35, 60, 180),
    "bone_marrow": _region("Bone and marrow", "rest", 0.050, 180, 2200, 5.0, 0.40, 48, 240),
    "rest": _region("Other tissues", "rest", 0.110, 300, 4200, 5.0, 0.45, 48, 180),
    "lung_left": _region("Left lung", "lung", 0.450, 245, 330, 35.0, 0.80, 16, 180, "pulmonary"),
    "lung_right": _region("Right lung", "lung", 0.550, 300, 400, 35.0, 0.80, 16, 180, "pulmonary"),
}

SYSTEMIC_REGIONS = tuple(key for key, value in REGIONS.items() if value.route == "systemic")
PULMONARY_REGIONS = tuple(key for key, value in REGIONS.items() if value.route == "pulmonary")

SYSTEMIC_FLOW_SUM = sum(REGIONS[key].flow_fraction for key in SYSTEMIC_REGIONS)
PULMONARY_FLOW_SUM = sum(REGIONS[key].flow_fraction for key in PULMONARY_REGIONS)
if not np.isclose(SYSTEMIC_FLOW_SUM, 1.0) or not np.isclose(PULMONARY_FLOW_SUM, 1.0):
    raise RuntimeError("Human regional blood-flow fractions must sum to one cardiac output")


CIRCULATION_VOLUMES_ML = {
    "arm_vein": 20.0,
    "right_heart": 100.0,
    "pulmonary_artery": 150.0,
    "pulmonary_vein": 150.0,
    "left_heart": 100.0,
    # Organ vascular spaces below already contain regional blood. Keeping the
    # central arterial/venous pools compact avoids counting that blood twice;
    # total central + regional vascular volume is approximately 5.06 L.
    "arterial": 500.0,
    "venous": 1000.0,
}

CIRCULATION_STATES = tuple(f"A_{name}" for name in CIRCULATION_VOLUMES_ML)
ROUTE_STATES = ("A_csf_lumbar", "A_csf_cranial", "A_im_depot", "A_airway_depot")
REGION_STATE_SUFFIXES = (
    "v", "isf", "bound", "ee", "le", "cy", "ncap", "nss", "nds",
    "epi", "mrna", "protein",
)
STATE_NAMES = list(CIRCULATION_STATES) + list(ROUTE_STATES)
for region_id in REGIONS:
    STATE_NAMES.extend(f"A_{region_id}_{suffix}" for suffix in REGION_STATE_SUFFIXES)
STATE_NAMES.extend(("Loss", "Dose_in"))
IDX = {name: index for index, name in enumerate(STATE_NAMES)}


def make_time_grid(max_days: float = 730.0) -> np.ndarray:
    """Multiscale grid resolving seconds, first pass, tissue PK, and persistence."""
    seconds = np.r_[0.0, np.geomspace(1.0 / 3600.0, 10.0 / 60.0, 52)]
    early_h = np.geomspace(10.5 / 60.0, 24.0, 82)
    days_h = np.geomspace(24.5, 30.0 * 24.0, 72)
    long_h = np.geomspace(31.0 * 24.0, max_days * 24.0, 80)
    return np.unique(np.r_[seconds, early_h, days_h, long_h])


def dose_input_rate(t_h: float, dose_vg: float, route: AdministrationRoute) -> float:
    duration = route.infusion_duration_h
    return dose_vg / duration if 0.0 <= t_h <= duration else 0.0


def solve_human_capsid(
    tropism: Mapping[str, float],
    dose_vg: float = DOSE_VG,
    t_eval: np.ndarray | None = None,
    administration: str = "iv",
):
    """Solve one capsid in the reference-adult regional model."""
    if administration not in ADMINISTRATION_ROUTES:
        raise ValueError(f"Unknown administration route: {administration}")
    administration_route = ADMINISTRATION_ROUTES[administration]
    times = make_time_grid() if t_eval is None else np.asarray(t_eval, dtype=float)
    y0 = np.zeros(len(STATE_NAMES), dtype=float)
    effective_cardiac_output = CARDIAC_OUTPUT_ML_H * EFFECTIVE_FLOW_SCALE
    k_blood_clear = np.log(2.0) / BLOOD_CAPSID_HALF_LIFE_H
    k_off = 0.05
    k_ee_le = 0.30
    k_rec = 0.05
    k_deg_ee = 0.02
    k_lys = 0.10
    k_nuc = 0.02
    k_uncoat_cyto = 0.005
    k_uncoat_nuc = 0.02
    k_deg_cyto = 0.01
    k_deg_ncap = 0.005
    k_ds = 0.01
    k_deg_ss = 0.02
    k_epi = 0.01
    k_deg_ds = 0.005
    k_mrna_deg = np.log(2.0) / 6.0
    k_protein_deg = np.log(2.0) / 48.0
    arm_transit_h = 18.0 / 3600.0

    def rhs(t_h: float, y: np.ndarray) -> np.ndarray:
        state = np.maximum(y, 0.0)
        d = np.zeros_like(state)
        input_rate = dose_input_rate(t_h, dose_vg, administration_route)

        def amount(name: str) -> float:
            return state[IDX[name]]

        def concentration(circulation: str) -> float:
            return amount(f"A_{circulation}") / CIRCULATION_VOLUMES_ML[circulation]

        c_rh = concentration("right_heart")
        c_pa = concentration("pulmonary_artery")
        c_pv = concentration("pulmonary_vein")
        c_lh = concentration("left_heart")
        c_art = concentration("arterial")
        c_ven = concentration("venous")
        arm_out = amount("A_arm_vein") / arm_transit_h

        d[IDX[administration_route.input_state]] += input_rate
        d[IDX["A_arm_vein"]] -= arm_out
        d[IDX["A_right_heart"]] += arm_out + effective_cardiac_output * (c_ven - c_rh)
        d[IDX["A_pulmonary_artery"]] += effective_cardiac_output * (c_rh - c_pa)
        d[IDX["A_left_heart"]] += effective_cardiac_output * (c_pv - c_lh)
        d[IDX["A_arterial"]] += effective_cardiac_output * (c_lh - c_art)

        pulmonary_return = 0.0
        venous_return = 0.0
        portal_gut = 0.0
        portal_spleen = 0.0

        # Lumbar IT dosing first occupies the caudal CSF, then moves rostrally.
        # Both CSF pools can drain to venous blood, preserving the peripheral
        # exposure observed after CSF administration.
        c_csf_lumbar = amount("A_csf_lumbar") / ROUTE_COMPARTMENT_VOLUMES_ML["csf_lumbar"]
        c_csf_cranial = amount("A_csf_cranial") / ROUTE_COMPARTMENT_VOLUMES_ML["csf_cranial"]
        if administration == "intrathecal":
            lumbar_to_cranial = CSF_ROSTRAL_FLOW_ML_H * c_csf_lumbar
            cranial_to_lumbar = 0.0
        elif administration in {"intracisternal", "intracerebroventricular"}:
            lumbar_to_cranial = 0.0
            cranial_to_lumbar = 0.60 * CSF_ROSTRAL_FLOW_ML_H * c_csf_cranial
        else:
            lumbar_to_cranial = 0.0
            cranial_to_lumbar = 0.0
        k_csf_absorb = np.log(2.0) / CSF_ABSORPTION_HALF_LIFE_H
        lumbar_absorption = 0.35 * k_csf_absorb * amount("A_csf_lumbar")
        cranial_absorption = 0.65 * k_csf_absorb * amount("A_csf_cranial")
        d[IDX["A_csf_lumbar"]] += cranial_to_lumbar - lumbar_to_cranial - lumbar_absorption
        d[IDX["A_csf_cranial"]] += lumbar_to_cranial - cranial_to_lumbar - cranial_absorption
        d[IDX["A_venous"]] += lumbar_absorption + cranial_absorption

        # IM dosing is represented as a local extracellular depot. Most vector
        # reaches injected muscle ISF; a smaller fraction escapes systemically.
        im_release = (np.log(2.0) / IM_DEPOT_ABSORPTION_HALF_LIFE_H) * amount("A_im_depot")
        d[IDX["A_im_depot"]] -= im_release
        d[IDX["A_muscle_injected_arm_isf"]] += IM_LOCAL_FRACTION * im_release
        d[IDX["A_arm_vein"]] += (1.0 - IM_LOCAL_FRACTION) * im_release

        # Airway delivery is a local pulmonary depot. The retained fraction is
        # split by lung ISF volume; the remainder reaches pulmonary venous blood.
        airway_release = (np.log(2.0) / AIRWAY_DEPOT_RELEASE_HALF_LIFE_H) * amount("A_airway_depot")
        d[IDX["A_airway_depot"]] -= airway_release
        lung_isf_total = REGIONS["lung_left"].isf_ml + REGIONS["lung_right"].isf_ml
        d[IDX["A_lung_left_isf"]] += AIRWAY_LOCAL_FRACTION * airway_release * REGIONS["lung_left"].isf_ml / lung_isf_total
        d[IDX["A_lung_right_isf"]] += AIRWAY_LOCAL_FRACTION * airway_release * REGIONS["lung_right"].isf_ml / lung_isf_total
        d[IDX["A_pulmonary_vein"]] += (1.0 - AIRWAY_LOCAL_FRACTION) * airway_release

        for region_id, region in REGIONS.items():
            v_name = f"A_{region_id}_v"
            isf_name = f"A_{region_id}_isf"
            bound_name = f"A_{region_id}_bound"
            ee_name = f"A_{region_id}_ee"
            le_name = f"A_{region_id}_le"
            cy_name = f"A_{region_id}_cy"
            ncap_name = f"A_{region_id}_ncap"
            nss_name = f"A_{region_id}_nss"
            nds_name = f"A_{region_id}_nds"
            epi_name = f"A_{region_id}_epi"
            mrna_name = f"A_{region_id}_mrna"
            protein_name = f"A_{region_id}_protein"
            c_v = amount(v_name) / region.vascular_ml
            c_isf = amount(isf_name) / region.isf_ml
            multiplier = max(float(tropism.get(region.parent, 1.0)), 0.01)
            ps = region.ps_ml_h * multiplier
            kp = region.kp * np.sqrt(multiplier)
            exchange = ps * (c_v - c_isf / max(kp, 1e-12))

            csf_exchange = 0.0
            csf_source_state = None
            if region_id in CSF_REGION_ACCESS:
                csf_pool, csf_ps, csf_kp = CSF_REGION_ACCESS[region_id]
                csf_ps *= CSF_ROUTE_ACCESS_MULTIPLIER.get(administration, {}).get(region_id, 1.0)
                csf_source_state = f"A_csf_{csf_pool}"
                c_csf = c_csf_lumbar if csf_pool == "lumbar" else c_csf_cranial
                csf_exchange = csf_ps * (c_csf - c_isf / csf_kp)
                if amount(csf_source_state) <= 0.0 and csf_exchange > 0.0:
                    csf_exchange = 0.0

            capsid_loss = np.log(2.0) / CAPSID_HALF_LIFE_H[region.parent]
            vascular_loss = 0.35 * capsid_loss * amount(v_name)
            isf_loss = 0.65 * capsid_loss * amount(isf_name)

            bound = amount(bound_name)
            receptor_capacity = RECEPTOR_DENSITY_VG_ML[region.parent] * region.isf_ml
            free_receptor = max(receptor_capacity - bound, 0.0)
            binding = BINDING_ON_ML_VG_H[region.parent] * c_isf * free_receptor - k_off * bound
            if bound <= 0.0 and binding < 0.0:
                binding = 0.0
            k_int = np.log(2.0) / region.internalization_half_life_h
            internalization = k_int * bound

            ee = amount(ee_name)
            le = amount(le_name)
            cy = amount(cy_name)
            ncap = amount(ncap_name)
            nss = amount(nss_name)
            nds = amount(nds_name)
            epi = amount(epi_name)
            k_escape = 0.008 if region.parent == "brain" else 0.006 if region.parent == "kidney" else 0.005
            epi_loss = (np.log(2.0) / (24.0 * region.episome_half_life_days)) * amount(epi_name)

            flow = effective_cardiac_output * region.flow_fraction
            if region.route == "pulmonary":
                flow_flux = flow * (c_pa - c_v)
                pulmonary_return += flow * c_v
            elif region_id == "liver":
                c_gut = amount("A_gut_v") / REGIONS["gut"].vascular_ml
                c_spleen = amount("A_spleen_v") / REGIONS["spleen"].vascular_ml
                hepatic_artery = flow * c_art
                portal_gut = effective_cardiac_output * REGIONS["gut"].flow_fraction * c_gut
                portal_spleen = effective_cardiac_output * REGIONS["spleen"].flow_fraction * c_spleen
                liver_out_flow = flow + effective_cardiac_output * (
                    REGIONS["gut"].flow_fraction + REGIONS["spleen"].flow_fraction
                )
                flow_flux = hepatic_artery + portal_gut + portal_spleen - liver_out_flow * c_v
                venous_return += liver_out_flow * c_v
            else:
                flow_flux = flow * (c_art - c_v)
                if region_id not in {"gut", "spleen"}:
                    venous_return += flow * c_v

            d[IDX[v_name]] += flow_flux - exchange - vascular_loss
            d[IDX[isf_name]] += exchange + csf_exchange - binding - isf_loss
            if csf_source_state is not None:
                d[IDX[csf_source_state]] -= csf_exchange
            d[IDX[bound_name]] += binding - internalization
            d[IDX[ee_name]] += internalization - (k_ee_le + k_rec + k_deg_ee) * ee
            d[IDX[le_name]] += k_ee_le * ee - (k_escape + k_lys) * le
            d[IDX[cy_name]] += k_escape * le - (k_nuc + k_uncoat_cyto + k_deg_cyto) * cy
            d[IDX[ncap_name]] += k_nuc * cy - (k_uncoat_nuc + k_deg_ncap) * ncap
            d[IDX[nss_name]] += k_uncoat_cyto * cy + k_uncoat_nuc * ncap - (k_ds + k_deg_ss) * nss
            d[IDX[nds_name]] += k_ds * nss - (k_epi + k_deg_ds) * nds
            d[IDX[epi_name]] += k_epi * nds - epi_loss

            h = 1.2
            ec50_tx = max(0.01 * receptor_capacity, 100.0)
            tx = 2.0 * epi**h / (ec50_tx**h + epi**h + 1e-30)
            d[IDX[mrna_name]] += min(tx, 2.0) - k_mrna_deg * amount(mrna_name)
            d[IDX[protein_name]] += 5.0 * amount(mrna_name) - k_protein_deg * amount(protein_name)

            intracellular_loss = (
                (k_rec + k_deg_ee) * ee
                + k_lys * le
                + k_deg_cyto * cy
                + k_deg_ncap * ncap
                + k_deg_ss * nss
                + k_deg_ds * nds
                + epi_loss
            )
            d[IDX["Loss"]] += vascular_loss + isf_loss + intracellular_loss

        d[IDX["A_pulmonary_vein"]] += pulmonary_return - effective_cardiac_output * c_pv
        d[IDX["A_venous"]] += venous_return - effective_cardiac_output * c_ven

        for circulation in CIRCULATION_VOLUMES_ML:
            if circulation == "arm_vein":
                continue
            cleared = k_blood_clear * amount(f"A_{circulation}")
            d[IDX[f"A_{circulation}"]] -= cleared
            d[IDX["Loss"]] += cleared

        d[IDX["Dose_in"]] = input_rate
        return d

    solution = solve_ivp(
        rhs,
        (float(times[0]), float(times[-1])),
        y0,
        t_eval=times,
        method="BDF",
        dense_output=True,
        rtol=2e-7,
        atol=1e-4,
        max_step=24.0,
    )
    if not solution.success:
        raise RuntimeError(f"Human spatial PBPK solve failed: {solution.message}")
    return solution


def vector_mass(solution) -> np.ndarray:
    total = np.zeros_like(solution.t)
    for name in CIRCULATION_STATES:
        total += solution.y[IDX[name]]
    for name in ROUTE_STATES:
        total += solution.y[IDX[name]]
    for region_id in REGIONS:
        for suffix in (
            "v", "isf", "bound", "ee", "le", "cy", "ncap", "nss", "nds", "epi",
        ):
            total += solution.y[IDX[f"A_{region_id}_{suffix}"]]
    return total


def mass_balance_error(solution) -> np.ndarray:
    delivered = np.maximum(solution.y[IDX["Dose_in"]], 1e-30)
    return (vector_mass(solution) + solution.y[IDX["Loss"]] - solution.y[IDX["Dose_in"]]) / delivered


def main() -> None:
    tropism = {parent: 1.0 for parent in {r.parent for r in REGIONS.values()}}
    print(f"states={len(STATE_NAMES)} regions={len(REGIONS)} routes={len(ADMINISTRATION_ROUTES)}")
    for route_id in ADMINISTRATION_ROUTES:
        solution = solve_human_capsid(tropism, administration=route_id)
        error = mass_balance_error(solution)
        meaningful = solution.y[IDX["Dose_in"]] > DOSE_VG * 1e-9
        accounted = (vector_mass(solution)[-1] + solution.y[IDX["Loss"], -1]) / solution.y[IDX["Dose_in"], -1]
        print(
            f"{route_id}: max_mass_balance_error={np.max(np.abs(error[meaningful])):.3e} "
            f"final_accounted_fraction={accounted:.6f}"
        )


if __name__ == "__main__":
    main()
