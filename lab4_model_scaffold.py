"""
MEGN 417/527 Vehicle Dynamics -- Lab 3
Numerical Longitudinal Model: Scaffold

Fill in every function marked TODO. The structure, the diagnostics, and the
validation harness are provided. The physics is not.

Run:  python lab3_model_scaffold.py
The validation suite runs automatically and reports PASS/FAIL for each test.
Tests will fail until you implement the model. That is expected.

Python 3.9+, numpy, scipy, matplotlib.
"""

import numpy as np
from scipy.integrate import solve_ivp

# numpy >=2.0 renamed trapz -> trapezoid; support both
_trapz = getattr(np, "trapezoid", None) or np.trapz
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
# Replace every value with your Specification Sheet v0.2 entry.
# The 'status' dict is not decoration. Keep it current. You will be asked
# which of your results depend on assumed rather than measured values.

P = {
    # --- mass ---
    "m_vehicle":  50.17,   # kg, unloaded, MEASURED in Lab 2
    "m_driver":    67.1,   # kg
    # --- geometry ---
    "L":        0.9017,   # m, wheelbase at your locked frame setting
    "a":          0.57,   # m, front axle to CG (from corner weights), loaded
    "h":        0.5735,   # m, combined CG height (Lab 2: tilt-measured kart
                          #      + anthropometric driver, composed) MAY NEED TO BE REVISED - Raistlin
    "h_aero":     0.5735,   # m, height at which drag acts
    "t_track":  0.7366,   # m, track width FRONT
    "t_track_r":0.50165,   # m, track width REAR
    "r_w":        0.127,   # m, rolling radius
    # --- rotating inertia ---
    "I_wheel_f":  0.008,   # kg*m^2, per front wheel
    "I_wheel_r":  0.016,   # kg*m^2, per rear wheel ESTIMATED WITH BOUNDING METHOD - Raistlin
    "I_motor":    0.0,   # kg*m^2, per motor rotor (0.0 if lumped into wheel)
    "N_gear":     1.0,    # motor:wheel ratio. 1.0 for direct-drive hub motors
    # --- resistances ---
    "C_rr":       0.012,   # - PROVISIONAL
    "C_d":        0.7,   # -
    "A_frontal":  0.25,   # m^2 - Provided by Colin Thomas
    "rho_air":    1,   # kg/m^3  -- check against Golden CO elevation!
    # --- tire / road ---
    "mu":         5,   # - (ASSUMED until Lab 4) Unsprung Mass Per Corner, PROVISIONAL
    "k_t_f":      1.316*1e6,   # N/m, front tire vertical stiffness
    "k_t_r":      8.225*1e5,   # N/m, rear tire vertical stiffness
    # --- powertrain ---
    "T_max":      96,   # N*m, total at the wheel
    "P_max":      4800,   # W, total. Published: nominal 2x450 W, peak 2x2400 W.
                          #    Run Test 5 with both and report both.
    "eta":        1.0,   # - driveline efficiency
    "drive_axle": "rear", # 'rear' | 'front' | 'awd'
    # --- constants ---
    "g":          9.81,
}

STATUS = {
    # parameter: 'measured' | 'derived' | 'bounded' | 'assumed' | 'published'
    # Fill this in. It drives your uncertainty discussion.
}


def total_mass(p):
    return p["m_vehicle"] + p["m_driver"]


# ---------------------------------------------------------------------------
# MODEL COMPONENTS  -- implement these
# ---------------------------------------------------------------------------

def effective_mass(p):
    """
    Return m_eff [kg].

    m_eff = m + (sum of wheel inertias + motor rotor inertia * N^2) / r_w^2

    Careful: m_eff applies to the INERTIAL term only. Do not use it for
    rolling resistance or grade.
    """
    m = total_mass(p)
    I_total = 2 * p["I_wheel_f"] + 2 * p["I_wheel_r"]
    I_total += 2 * p["I_motor"] * p["N_gear"] ** 2
    return m + I_total / p["r_w"] ** 2
    raise NotImplementedError


def inertia_fraction(p):
    """
    Return gamma = (m_eff - m) / m = sum(I) / (m r_w^2) [-].

    Course convention: gamma is the INCREMENT, so m_eff = m * (1 + gamma).
    Expect a few percent. Zero inertia gives gamma = 0 exactly.
    """
    m = total_mass(p)
    return (effective_mass(p) - m) / m
    raise NotImplementedError


def aero_drag(v, p):
    """Return aerodynamic drag [N], positive opposing motion."""
    # TODO
    return 0.5 * p["rho_air"] * p["C_d"] * p["A_frontal"] * v ** 2
    raise NotImplementedError


def rolling_resistance(v, p, grade=0.0):
    """
    Return rolling resistance [N], positive opposing motion.
    grade is the road angle in radians.
    Use ACTUAL mass here, not effective mass.
    """
    # TODO
    m = total_mass(p)
    return p["C_rr"] * m * p["g"] * np.cos(grade)
    raise NotImplementedError


def grade_force(p, grade=0.0):
    """Return grade resistance [N], positive opposing motion on an incline."""
    # TODO
    m = total_mass(p)
    return m * p["g"] * np.sin(grade)
    raise NotImplementedError


def rear_axle_load(a_x, v, p):
    """
    Return dynamic rear axle normal force [N].

        F_zr = (a/L) m g + (h/L) m a_x + (h_aero/L) F_aero

    Watch the sign of the a_x term. Forward acceleration must INCREASE
    rear load on this vehicle. If it decreases, you have the axle
    positions reversed relative to the direction of travel.
    """
    # TODO
    m = total_mass(p)
    F_aero = aero_drag(v, p)
    F_zr = (
        (p["a"] / p["L"]) * m * p["g"]
        + (p["h"] / p["L"]) * m * a_x
        + (p["h_aero"] / p["L"]) * F_aero
    )

    return F_zr
    raise NotImplementedError


def traction_limit(a_x, v, p):
    """Return max transmissible drive force [N] = mu * F_z on the driven axle."""
    # TODO
    rear = rear_axle_load(a_x, v, p)

    front = total_mass(p) * p["g"] - rear

    if p["drive_axle"] == "rear":
        normal = rear

    elif p["drive_axle"] == "front":
        normal = front

    elif p["drive_axle"] == "awd":
        normal = front + rear

    else:
        raise ValueError("Invalid drive_axle")

    return p["mu"] * max(normal, 0.0)
    raise NotImplementedError


def drive_force(v, a_x, p, throttle=1.0):
    """
    Return drive force [N] as the minimum of three ceilings:
        torque limit   : T_max * N * eta / r_w
        power limit    : P_max * eta / v     (guard v -> 0)
        traction limit : mu * F_zr(a_x, v)
    Scale the powertrain ceilings by throttle, not the traction ceiling.
    """
    # TODO
    torque_force = (
    throttle * p["T_max"] * p["N_gear"] * p["eta"] / p["r_w"]
    )

    power_force = (
        throttle * p["P_max"] * p["eta"] / max(v, 0.000001)
    )

    tire_force = traction_limit(a_x, v, p)

    return min(torque_force, power_force, tire_force)
    raise NotImplementedError


def binding_limit(v, a_x, p, throttle=1.0):
    """Return 'torque' | 'power' | 'traction'. Used for diagnostics."""
    # TODO
    ceilings = {
        "torque": throttle * p["T_max"] * p["N_gear"] * p["eta"] / p["r_w"],

        "power": throttle * p["P_max"] * p["eta"] / max(v, 0.000001),

        "traction": traction_limit(a_x, v, p),
    }

    return min(ceilings, key=ceilings.get)
    raise NotImplementedError


def accel(v, p, throttle=1.0, grade=0.0):
    """
    Return dv/dt [m/s^2].

    F_drive depends on a_x through load transfer, and a_x depends on
    F_drive. Resolve this properly:
      - closed form where drag allows, or
      - fixed-point iteration to a stated tolerance.

    Do NOT lag a_x by one timestep.

    Raise or flag if (m_eff - mu*m*h/L) <= 0: the front axle would lift
    before the tires slip, and the closed form is invalid. Liftoff can
    also bind earlier: check a_x against g*b/h.
    """
    # TODO
    m = total_mass(p)

    meff = effective_mass(p)

    resist = (
        aero_drag(v, p)
        + rolling_resistance(v, p, grade)
        + grade_force(p, grade)
    )

    force_powertrain = min(
        throttle * p["T_max"] * p["N_gear"] * p["eta"] / p["r_w"],
        throttle * p["P_max"] * p["eta"] / max(v, 0.000001),
    )

    a_candidate = (force_powertrain - resist) / meff

    if force_powertrain <= traction_limit(a_candidate, v, p):
        a_x = a_candidate

    else:
        drive_axle = p["drive_axle"]

        F_at_zero = traction_limit(0.0, v, p)

        if drive_axle == "rear":
            slope = p["mu"] * m * p["h"] / p["L"]

        elif drive_axle == "front":
            slope = -p["mu"] * m * p["h"] / p["L"]

        else:
            slope = 0.0

        denominator = meff - slope

        if denominator <= 0:
            raise ValueError(
                "Traction-limited calculation invalid: possible front liftoff"
            )

        a_x = (F_at_zero - resist) / denominator

    b = p["L"] - p["a"]

    if p["h"] > 0 and a_x >= p["g"] * b / p["h"]:
        raise ValueError("Front axle liftoff predicted")

    if total_mass(p) * p["g"] - rear_axle_load(a_x, v, p) <= 0:
        raise ValueError("Front axle normal load is nonpositive")

    return a_x
    raise NotImplementedError


# ---------------------------------------------------------------------------
# INTEGRATION
# ---------------------------------------------------------------------------

def simulate(p, v0=0.0, t_end=30.0, throttle=1.0, grade=0.0, rtol=1e-8, atol=1e-10):
    """
    Integrate the model. Returns a dict with t, v, x, a, and diagnostics.
    """
    def rhs(t, y):
        v = max(y[0], 0.0)
        return [accel(v, p, throttle, grade), v]

    sol = solve_ivp(rhs, [0, t_end], [v0, 0.0],
                    dense_output=True, rtol=rtol, atol=atol, max_step=0.05)

    t = np.linspace(0, t_end, 2000)
    v, x = sol.sol(t)
    v = np.maximum(v, 0.0)

    a = np.array([accel(vi, p, throttle, grade) for vi in v])
    Fzr = np.array([rear_axle_load(ai, vi, p) for ai, vi in zip(a, v)])
    Fzf = total_mass(p) * p["g"] - Fzr
    limit = [binding_limit(vi, ai, p, throttle) for vi, ai in zip(v, a)]
    power = np.array([drive_force(vi, ai, p, throttle) * vi
                      for vi, ai in zip(v, a)])

    return {"t": t, "v": v, "x": x, "a": a,
            "F_zr": Fzr, "F_zf": Fzf, "limit": limit, "power": power,
            "sol": sol}


# ---------------------------------------------------------------------------
# VALIDATION SUITE  -- do not edit; your model must satisfy these
# ---------------------------------------------------------------------------

def _report(name, passed, detail=""):
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}"
          + (f"   {detail}" if detail else ""))
    return passed


def validate(p):
    print("\nVALIDATION SUITE")
    print("-" * 70)
    results = []

    # Test 1: zero inertia collapses to gamma = 0 and m_eff = m
    try:
        q = dict(p)
        q["I_wheel_f"] = q["I_wheel_r"] = q["I_motor"] = 0.0
        g1 = inertia_fraction(q)
        ok = abs(g1) < 1e-12 and abs(effective_mass(q) - total_mass(q)) < 1e-9
        results.append(_report("1. gamma -> 0 with zero inertia", ok,
                              f"gamma = {g1:.12f}"))
    except Exception as e:
        results.append(_report("1. gamma -> 0 with zero inertia", False, str(e)))

    # Test 2: a = F/m_eff with all resistances off and load transfer disabled.
    # Setting h = h_aero = 0 removes the a_x dependence of F_zr, so the drive
    # force no longer depends on a_x and the comparison is well posed whichever
    # ceiling binds. Do not instead set mu very large: that trips the wheelie
    # guard, because the guard tests gamma - mu*h/L.
    try:
        q = dict(p)
        q["C_rr"] = 0.0
        q["C_d"] = 0.0
        q["h"] = 0.0
        q["h_aero"] = 0.0
        m_eff = effective_mass(q)
        F = drive_force(5.0, 0.0, q)
        a_expect = F / m_eff
        a_actual = accel(5.0, q, throttle=1.0, grade=0.0)
        ok = abs(a_actual - a_expect) / max(abs(a_expect), 1e-9) < 1e-6
        results.append(_report("2. a = F/m_eff, no resistances", ok,
                              f"expected {a_expect:.6f}, got {a_actual:.6f}"))
    except Exception as e:
        results.append(_report("2. a = F/m_eff, no resistances", False, str(e)))

    # Test 3: vertical equilibrium
    try:
        r = simulate(p, t_end=10.0)
        W = total_mass(p) * p["g"]
        err = np.max(np.abs(r["F_zf"] + r["F_zr"] - W))
        results.append(_report("3. F_zf + F_zr = mg", err < 1e-6,
                              f"max error {err:.3e} N"))
    except Exception as e:
        results.append(_report("3. F_zf + F_zr = mg", False, str(e)))

    # Test 4: energy audit
    try:
        r = simulate(p, t_end=15.0)
        t, v, a = r["t"], r["v"], r["a"]
        m_eff = effective_mass(p)
        Fd = np.array([drive_force(vi, ai, p) for vi, ai in zip(v, a)])
        Fr = np.array([rolling_resistance(vi, p) for vi in v])
        Fa = np.array([aero_drag(vi, p) for vi in v])
        W_in = _trapz(Fd * v, t)
        W_loss = _trapz((Fr + Fa) * v, t)
        dKE = 0.5 * m_eff * (v[-1] ** 2 - v[0] ** 2)
        rel = abs(W_in - dKE - W_loss) / max(abs(W_in), 1e-9)
        results.append(_report("4. energy audit closes", rel < 1e-3,
                              f"residual {rel*100:.4f}% of input work"))
    except Exception as e:
        results.append(_report("4. energy audit closes", False, str(e)))

    # Test 5: terminal velocity vs published top speed
    try:
        r = simulate(p, t_end=90.0)
        v_term = r["v"][-1]
        print(f"  [INFO] 5. terminal velocity {v_term:.2f} m/s "
              f"= {v_term*3.6:.1f} km/h  -- compare to published. "
              f"Do NOT tune to match.")
    except Exception as e:
        print(f"  [INFO] 5. terminal velocity: {e}")

    # Test 6: solver independence
    try:
        r1 = simulate(p, t_end=15.0, rtol=1e-8)
        r2 = simulate(p, t_end=15.0, rtol=1e-10)
        d = np.max(np.abs(r1["v"] - r2["v"])) / max(np.max(r1["v"]), 1e-9)
        results.append(_report("6. solver independence", d < 1e-3,
                              f"max relative change {d*100:.5f}%"))
    except Exception as e:
        results.append(_report("6. solver independence", False, str(e)))

    print("-" * 70)
    print(f"  {sum(results)}/{len(results)} tests passed\n")
    return results


# ---------------------------------------------------------------------------
# REQUIRED RUNS
# ---------------------------------------------------------------------------

def required_runs(p):
    """Produce the four runs required by the handout."""
    # Run 1: full-throttle launch
    # Run 2: coastdown from top speed (throttle=0, v0=v_top)
    # Run 3: constant-speed power demand, 5..45 km/h
    # Run 4: launch on the manufacturer's max rated grade
    # TODO: implement and plot. Label every axis with units.
    
    #Run 1
    # RUN 1

    try:
        launch = simulate(p, t_end=30.0)

        plt.figure()
        plt.plot(launch["t"], launch["v"] * 3.6)
        plt.xlabel("Time (s)")
        plt.ylabel("Speed (km/h)")
        plt.title("Run 1: Full-throttle launch")
        plt.grid(True)

    except ValueError as e:
        print(f"Run 1 failed: {e}")


    # RUN 2

    v_top = None

    try:
        top_run = simulate(p, t_end=90.0)
        v_top = top_run["v"][-1]

        coast = simulate(p, v0=v_top, t_end=30.0, throttle=0.0)

        plt.figure()
        plt.plot(coast["t"], coast["v"] * 3.6)
        plt.xlabel("Time (s)")
        plt.ylabel("Speed (km/h)")
        plt.title("Run 2: Coastdown")
        plt.grid(True)

    except ValueError as e:
        print(f"Run 2 failed: {e}")

    # RUN 3

    speeds_kmh = np.linspace(5.0, 45.0, 100)

    speeds_ms = speeds_kmh / 3.6

    forces = np.array([
        aero_drag(vi, p) + rolling_resistance(vi, p)
        for vi in speeds_ms
    ])

    wheel_power = forces * speeds_ms

    plt.figure()

    plt.plot(speeds_kmh, wheel_power)

    plt.xlabel("Speed (km/h)")
    plt.ylabel("Required wheel power (W)")
    plt.title("Run 3: Constant-speed power demand")
    plt.grid(True)


    # RUN 4

    rated_grade_percent = float(
        input("Manufacturer maximum rated grade (%): ")
    )

    grade_angle = np.arctan(rated_grade_percent / 100.0)

    try:
        incline = simulate(p, t_end=30.0, grade=grade_angle)

        plt.figure()
        plt.plot(incline["t"], incline["v"] * 3.6)
        plt.xlabel("Time (s)")
        plt.ylabel("Speed (km/h)")
        plt.title(f"Run 4: Launch on {rated_grade_percent:g}% grade")
        plt.grid(True)

    except ValueError as e:
        print(f"Run 4 failed: {e}")


    # NOMINAL VS PEAK POWER COMPARISON

    nominal_p = dict(p)

    nominal_p["P_max"] = 900.0

    try:
        nominal_top = simulate(nominal_p, t_end=90.0)["v"][-1]
        print(f"Nominal-power 90 s speed: {nominal_top * 3.6:.2f} km/h")

    except ValueError as e:
        print(f"Nominal-power simulation failed: {e}")

    if v_top is not None:
        print(f"Peak-power 90 s speed: {v_top * 3.6:.2f} km/h")
    else:
        print("Peak-power simulation failed.")


    # DISPLAY FIGURES

    plt.show()

    return None
    raise NotImplementedError


if __name__ == "__main__":
    missing = [k for k, v in P.items() if v is None]
    if missing:
        print("\nFill in these parameters before running:")
        for k in missing:
            print(f"   P['{k}']")
        print()
    else:
        validate(P)
        required_runs(P)
