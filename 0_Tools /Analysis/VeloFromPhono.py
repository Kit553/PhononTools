#!/usr/bin/env python3


from pathlib import Path
import numpy as np
import phonopy

# PATH SETTINGS
PHONOPY_DIR = Path(r"path/to/PhonoPy/results")
OUTPUT_DIR = Path(r"path/to/output/directory")

YAML_FILE = "phonopy_disp.yaml" 
FORCE_SETS_FILE = "FORCE_SETS"

# NAC
USE_NAC = True
BORN_FILE = "BORN"

MESH = [20, 20, 20]    # Adjust this to converge the group velocities

# near-Gamma shell in reduced coordinates
Q_MIN = 1e-6
Q_MAX = 0.18

# Phonopy group velocities for VASP are usually THz Å.
# 1 THz Å = 100 m/s
GV_TO_M_S = 100.0


def harmonic_sound_velocity(v_t, v_l):
    return (3.0 / (v_l**-3 + 2.0 * v_t**-3)) ** (1.0 / 3.0)


def wrap_to_gamma(q):
    return q - np.round(q)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    yaml_path = PHONOPY_DIR / YAML_FILE
    force_sets_path = PHONOPY_DIR / FORCE_SETS_FILE
    born_path = PHONOPY_DIR / BORN_FILE

    print(f"Loading: {yaml_path}")
    print(f"Forces:  {force_sets_path}")

    load_kwargs = {
        "phonopy_yaml": str(yaml_path),
        "force_sets_filename": str(force_sets_path),
        "log_level": 1,
    }

    if USE_NAC:
        if born_path.exists():
            print(f"NAC:     using {born_path}")
            load_kwargs["born_filename"] = str(born_path)
            load_kwargs["is_nac"] = True
        else:
            print(f"NAC:     requested, but no BORN file found at {born_path}")
            print("         continuing without NAC")
            load_kwargs["is_nac"] = False
    else:
        print("NAC:     disabled")
        load_kwargs["is_nac"] = False

    phonon = phonopy.load(**load_kwargs)

    if phonon.force_constants is None:
        phonon.produce_force_constants()

    phonon.run_mesh(
        MESH,
        with_group_velocities=True,
        is_mesh_symmetry=False,
    )

    md = phonon.get_mesh_dict()

    qpoints = np.asarray(md["qpoints"])
    freqs = np.asarray(md["frequencies"])
    gv = np.asarray(md["group_velocities"]) * GV_TO_M_S

    q_wrapped = wrap_to_gamma(qpoints)
    q_dist = np.linalg.norm(q_wrapped, axis=1)

    mask = (q_dist >= Q_MIN) & (q_dist <= Q_MAX)

    if mask.sum() == 0:
        raise RuntimeError("No near-Gamma q-points found. Increase Q_MAX.")

    # first three branches = acoustic
    speeds = np.linalg.norm(gv[mask, :3, :], axis=2)
    freqs_ac = freqs[mask, :3]
    q_sel = q_wrapped[mask]

    v_l = []
    v_t = []
    rows = []

    for q, f3, s3 in zip(q_sel, freqs_ac, speeds):
        order = np.argsort(s3)

        vt1 = s3[order[0]]
        vt2 = s3[order[1]]
        vl = s3[order[2]]

        v_t.extend([vt1, vt2])
        v_l.append(vl)

        rows.append([
            q[0], q[1], q[2],
            f3[0], f3[1], f3[2],
            s3[0], s3[1], s3[2],
            vt1, vt2, vl,
        ])

    v_l = np.asarray(v_l)
    v_t = np.asarray(v_t)

    vL_mean = np.mean(v_l)
    vT_mean = np.mean(v_t)
    vS = harmonic_sound_velocity(vT_mean, vL_mean)

    out = OUTPUT_DIR / "near_gamma_velocities.csv"

    header = (
        "q1,q2,q3,"
        "freq1_THz,freq2_THz,freq3_THz,"
        "speed1_m_s,speed2_m_s,speed3_m_s,"
        "vT1_sorted_m_s,vT2_sorted_m_s,vL_sorted_m_s"
    )

    np.savetxt(
        out,
        np.asarray(rows),
        delimiter=",",
        header=header,
        comments="",
    )

    print()
    print("==========================================")
    print("Near-Gamma acoustic velocity estimate")
    print("==========================================")
    print(f"Selected q-points: {mask.sum()}")
    print(f"v_L = {vL_mean:.1f} ± {np.std(v_l):.1f} m/s")
    print(f"v_T = {vT_mean:.1f} ± {np.std(v_t):.1f} m/s")
    print(f"v_s = {vS:.1f} m/s")
    print()
    print("Paste into fit script:")
    print(f"v_longitudinal_m_s={vL_mean:.1f},")
    print(f"v_transverse_m_s={vT_mean:.1f},")
    print()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
