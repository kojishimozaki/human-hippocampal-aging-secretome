#!/usr/bin/env python
"""receiver-map STEP 00 — freeze the receptor universe (manifest-authoritative) + verify the SHA pin.

Pre-registration: docs/RECEIVER_MAP_PLAN.md, tag prereg/receiver-map-v1.
env: bio.  Reads only committed/SHA-pinned inputs; writes results/receiver/.

The 22 age-UP ligands and the 76-receptor universe are reused VERBATIM from the frozen CONSEQUENCE
manifest (results/resilience/frozen_receiver_set_manifest.json). The 76 = frozen_receptors(56) U
excluded_receptors(20) [pre-detectability universe]. We do NOT re-derive or re-test the ligand set.
The ligand->receptor edges (for n_ligands annotation + LOO readiness) are extracted from the
SHA-pinned NicheNet rds in the companion R step 00b, cross-checked against this manifest 76.
"""
import os, sys, json, hashlib, pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
os.chdir(PROJ)
OUT = "results/receiver"; os.makedirs(OUT, exist_ok=True)

MANIFEST = "results/resilience/frozen_receiver_set_manifest.json"
RDS = "refs/nichenet/lr_network_human_21122021.rds"
RDS_SHA_PIN = "47c971d2fbba4ecd0ba7485d1846a74054432a0d1a97ea3e6a79ae27d0da8094"

m = json.load(open(MANIFEST))

ligands = sorted(m["ligands"])
assert len(ligands) == 22, f"expected 22 ligands, got {len(ligands)}"

universe76 = sorted(set(m["frozen_receptors"]) | set(m["excluded_receptors"]))
assert len(universe76) == 76, f"expected 76 receptors (56 U 20), got {len(universe76)}"
assert m["n_up_ligands"] == 22 and m["n_receptors_prefilter"] == 76

# SHA pin verification (rds is the ligand->receptor source of truth)
sha = hashlib.sha256(open(RDS, "rb").read()).hexdigest()
sha_ok = (sha == RDS_SHA_PIN)
assert sha_ok, f"rds SHA mismatch: {sha} != {RDS_SHA_PIN} (STOP: pre-reg integrity)"

pd.DataFrame({"ligand": ligands}).to_csv(f"{OUT}/ligands_22.csv", index=False)
pd.DataFrame({"receptor": universe76}).to_csv(f"{OUT}/_universe76_from_manifest.csv", index=False)

prov = {
    "step": "00_freeze_universe",
    "prereg_tag": "prereg/receiver-map-v1",
    "manifest": MANIFEST,
    "manifest_freeze_commit": m.get("git_commit"),
    "n_ligands": len(ligands),
    "n_receptors_universe": len(universe76),
    "rds": RDS, "rds_sha256": sha, "rds_sha_pin_ok": sha_ok,
}
json.dump(prov, open(f"{OUT}/_step00_provenance.json", "w"), indent=2)

print("STEP 00 freeze-universe")
print(f"  ligands (frozen, reused verbatim): {len(ligands)}")
print(f"  receptor universe 76 = frozen_receptors(56) U excluded_receptors(20): {len(universe76)}")
print(f"  rds SHA pin OK: {sha_ok}  ({sha[:16]}...)")
print(f"  wrote: {OUT}/ligands_22.csv, {OUT}/_universe76_from_manifest.csv, {OUT}/_step00_provenance.json")
