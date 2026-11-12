#!/usr/bin/env python3
"""Tabulate Stage 2 solver scans (data/stage2/scan_*.json): selection table at 1 amu and per-mass table. Usage: stage2_table.py <json> [<json> ...]"""
import json, sys
R = [r for f in sys.argv[1:] for r in json.load(open(f))]
ok = lambda r: r["margin_kT"] >= 2 and r["ratio_q_over_c"] >= 5 and 4 <= r["Bcl_PMF_kT"] <= 6
okv = lambda r: r["margin_kT"] >= 2 and r["ratio_q_over_c"] >= 5 and 4 <= r["V0_kT"] <= 6
print("| V0 | sigma | m | E0 | W(0) | margin (kT) | V0/kT | PMF barrier (kT) | rho_q/rho_c at 0 | rho_q(0)/rho_q,max | split (kcal) | f(V>E) gs / thermal | pass (PMF) | pass (V0/kT) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in sorted(R, key=lambda r: (r["m"] != 1.0, r["V0"], r["sigma"], -r["m"])):
    print("| %.1f | %.2f | %.1f | %.3f | %.3f | %.2f | %.2f | %.2f | %.2f | %.3f | %s | %.3f / %.3f | %s | %s |" % (
        r["V0"], r["sigma"], r["m"], r["E0"], r["W0"], r["margin_kT"], r["V0_kT"], r["Bcl_PMF_kT"], r["ratio_q_over_c"], r["rho_q0_over_qmax"],
        "%.3f" % r["split_kcal"] if r["split_kcal"] is not None else "-", r["f_ground"], r["f_thermal"],
        ("yes" if ok(r) else "no") if r["m"] == 1.0 else "", ("yes" if okv(r) else "no") if r["m"] == 1.0 else ""))
