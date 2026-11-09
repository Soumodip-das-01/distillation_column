"""Preliminary multicomponent debutanizer design.

This is an educational/preliminary design model based on assumed refinery data.
It uses Fenske, Underwood and Gilliland correlations. Replace the relative
volatilities and product recoveries with simulator or plant data for detailed
design.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Component:
    name: str
    z: float                 # feed mole fraction
    alpha: float             # relative volatility, relative to n-hexane
    mw: float                # kg/kmol
    distillate_recovery: float


COMPONENTS = (
    Component("n-Butane", 0.10, 3.50, 58.12, 1.0000),
    Component("n-Pentane", 0.20, 1.80, 72.15, 0.9500),
    Component("n-Hexane", 0.30, 1.00, 86.18, 0.0200),
    Component("n-Heptane", 0.25, 0.55, 100.21, 0.0010),
    Component("n-Octane", 0.15, 0.28, 114.23, 0.0005),
)

FEED_KMOL_H = 100.0
Q_FEED = 1.0                    # saturated liquid
TRAY_EFFICIENCY = 0.70
REFLUX_MULTIPLIER = 1.50
OPERATING_PRESSURE_KPA = 500.0
VAPOR_TEMPERATURE_K = 473.15
ALLOWABLE_VAPOR_VELOCITY = 1.30
TRAY_SPACING_M = 0.45
LATENT_HEAT_KJ_PER_KMOL = 30_000.0


def product_flows(feed_kmol_h: float):
    """Return component flows and total distillate/bottoms flows."""
    d = [feed_kmol_h * c.z * c.distillate_recovery for c in COMPONENTS]
    b = [feed_kmol_h * c.z - d_i for c, d_i in zip(COMPONENTS, d)]
    return d, b, sum(d), sum(b)


def normalized(values):
    total = sum(values)
    return [value / total for value in values]


def fenske_minimum_stages(x_d, x_b, light_key=1, heavy_key=2):
    """Minimum theoretical stages at total reflux."""
    lk = COMPONENTS[light_key]
    hk = COMPONENTS[heavy_key]
    separation_factor = (x_d[light_key] / x_d[heavy_key]) * (
        x_b[heavy_key] / x_b[light_key]
    )
    return math.log(separation_factor) / math.log(lk.alpha / hk.alpha)


def underwood_root(z, q=1.0, light_key=1, heavy_key=2):
    """Find the Underwood root between alpha_HK and alpha_LK by bisection."""
    alpha_lk = COMPONENTS[light_key].alpha
    alpha_hk = COMPONENTS[heavy_key].alpha

    def equation(theta):
        return sum(q * zi / (c.alpha - theta) for zi, c in zip(z, COMPONENTS)) - 1.0

    # The interval is deliberately kept away from the singular endpoints.
    lo = alpha_hk + 1e-8
    hi = alpha_lk - 1e-8
    if equation(lo) * equation(hi) > 0:
        raise ValueError("Could not bracket the Underwood root.")

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if abs(equation(mid)) < 1e-10:
            return mid
        if equation(lo) * equation(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def underwood_minimum_reflux(z, x_d, q=1.0):
    theta = underwood_root(z, q)
    r_plus_one = sum(
        xd * c.alpha / (c.alpha - theta)
        for xd, c in zip(x_d, COMPONENTS)
    )
    return theta, r_plus_one - 1.0


def gilliland_stages(n_min, reflux, reflux_min):
    """Estimate theoretical stages using the Eduljee form of Gilliland."""
    x = (reflux - reflux_min) / (reflux + 1.0)
    exponent = ((1.0 + 54.4 * x) / (11.0 + 117.2 * x)) * (
        (x - 1.0) / math.sqrt(x)
    )
    y = 1.0 - math.exp(exponent)
    return (n_min + y) / (1.0 - y), x, y


def design():
    z = [c.z for c in COMPONENTS]
    d_flow, b_flow, d_total, b_total = product_flows(FEED_KMOL_H)
    x_d = normalized(d_flow)
    x_b = normalized(b_flow)

    n_min = fenske_minimum_stages(x_d, x_b)
    theta, r_min = underwood_minimum_reflux(z, x_d, Q_FEED)
    reflux = REFLUX_MULTIPLIER * r_min
    n_theoretical, _, _ = gilliland_stages(n_min, reflux, r_min)
    n_actual = math.ceil(n_theoretical / TRAY_EFFICIENCY)

    liquid_reflux = reflux * d_total
    vapor_flow = (reflux + 1.0) * d_total
    vapor_mw = sum(x * c.mw for x, c in zip(x_d, COMPONENTS))
    vapor_mass_flow = vapor_flow * vapor_mw

    gas_constant = 8.314  # kPa m3/(kmol K)
    vapor_volume_m3_s = (
        vapor_flow * gas_constant * VAPOR_TEMPERATURE_K
        / OPERATING_PRESSURE_KPA / 3600.0
    )
    area_m2 = vapor_volume_m3_s / ALLOWABLE_VAPOR_VELOCITY
    diameter_m = math.sqrt(4.0 * area_m2 / math.pi)
    selected_diameter_m = math.ceil(diameter_m / 0.05) * 0.05
    column_height_m = n_actual * TRAY_SPACING_M + 3.5

    condenser_mw = vapor_flow * LATENT_HEAT_KJ_PER_KMOL / 3600.0 / 1000.0
    reboiler_mw = 1.10 * condenser_mw

    return locals()


def print_report(r):
    print("PRELIMINARY MULTICOMPONENT DEBUTANIZER DESIGN")
    print("=" * 52)
    print(f"Feed rate                 : {FEED_KMOL_H:,.2f} kmol/h")
    print(f"Distillate rate           : {r['d_total']:,.2f} kmol/h")
    print(f"Bottoms rate              : {r['b_total']:,.2f} kmol/h")
    print(f"Minimum stages (Fenske)   : {r['n_min']:.2f}")
    print(f"Underwood root            : {r['theta']:.4f}")
    print(f"Minimum reflux ratio      : {r['r_min']:.2f}")
    print(f"Selected reflux ratio     : {r['reflux']:.2f}")
    print(f"Theoretical stages        : {r['n_theoretical']:.2f}")
    print(f"Actual trays, 70% efficient: {r['n_actual']}")
    print(f"Reflux liquid flow        : {r['liquid_reflux']:.2f} kmol/h")
    print(f"Vapor flow                : {r['vapor_flow']:.2f} kmol/h")
    print(f"Estimated vapor flow      : {r['vapor_volume_m3_s']:.3f} m3/s")
    print(f"Calculated diameter       : {r['diameter_m']:.3f} m")
    print(f"Selected diameter         : {r['selected_diameter_m']:.2f} m")
    print(f"Estimated column height   : {r['column_height_m']:.2f} m")
    print(f"Condenser duty            : {r['condenser_mw']:.2f} MW")
    print(f"Reboiler duty             : {r['reboiler_mw']:.2f} MW")

    print("\nPRODUCT COMPOSITIONS (mole fraction)")
    print(f"{'Component':<12}{'Feed':>10}{'Distillate':>14}{'Bottoms':>12}")
    for c, z, xd, xb in zip(COMPONENTS, [c.z for c in COMPONENTS], r['x_d'], r['x_b']):
        print(f"{c.name:<12}{z:>10.4f}{xd:>14.4f}{xb:>12.4f}")


if __name__ == "__main__":
    print_report(design())
