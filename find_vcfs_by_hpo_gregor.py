#!/usr/bin/env python3
"""
Usage: python3 find_vcfs_by_hpo_gregor.py HP:0001250
"""

import sys
import firecloud.api as fapi

# workspaces we want to check
WORKSPACES = [
    ("gregor-ga4k", "GREGOR_GA4K1", "GA4K"),
    ("gregor-dcc", "GREGOR_COMBINED_CONSORTIUM_U12", "DCC"),
]

def fetch(ns, ws, table):
    r = fapi.get_entities(ns, ws, table)
    return r.json() if r.status_code == 200 else []

def vcf_type(path):
    p = path.lower()
    if ".sv." in p: return "SV"
    if ".cnv." in p: return "CNV"
    if "hard-filtered" in p: return "SNV"
    if "clinical_exome" in p: return "clinical"
    if "deepvariant" in p: return "DeepVariant"
    return "other"

def assay_type(path):
    if "/WGS/" in path or "_WGS" in path: return "WGS"
    if "/ES/" in path or "Exome" in path: return "Exome"
    if "/MGI/" in path or "_MGI" in path: return "MGI"
    if "/WGBS/" in path or "_WGBS" in path: return "WGBS"
    if "long-read" in path or "PacBio" in path: return "PacBio"
    return "unknown"


def find_vcfs(HPO):
    results = []
    
    for ns, ws_name, label in WORKSPACES:
        
        # who has this HPO?
        phenos = fetch(ns, ws_name, "phenotype")
        affected = set()
        for p in phenos:
            a = p.get("attributes", {})
            if a.get("term_id") == HPO:
                affected.add(a.get("participant_id"))
        
        if not affected:
            continue
        
        # grab participant info for family structure
        partis = fetch(ns, ws_name, "participant")
        parti_map = {}
        for p in partis:
            pid = p["name"]
            a = p.get("attributes", {})
            parti_map[pid] = {
                "fam": a.get("family_id", "?"),
                "mom": a.get("maternal_id") or "0",
                "dad": a.get("paternal_id") or "0",
            }
        
        # find trios - proband has HPO and both parents in dataset
        trios = []
        for pid in affected:
            if pid not in parti_map:
                continue
            info = parti_map[pid]
            mom, dad = info["mom"], info["dad"]
            
            trio = {
                "fam": info["fam"],
                "proband": pid,
                "mom": mom if mom != "0" and mom in parti_map else None,
                "dad": dad if dad != "0" and dad in parti_map else None,
            }
            trios.append(trio)
        
        # get all the vcf paths
        variants = fetch(ns, ws_name, "called_variants_dna_short_read")
        
        vcfs = {}  # set_id > list of vcf paths
        for v in variants:
            a = v.get("attributes", {})
            set_id = a.get("aligned_dna_short_read_set_id", "")
            vcf = a.get("called_variants_dna_file", "")
            if not vcf:
                continue
            if set_id not in vcfs:
                vcfs[set_id] = []
            vcfs[set_id].append(vcf)
        
        # match vcfs to trio members
        for trio in trios:
            members = [
                ("proband", trio["proband"]),
                ("mother", trio["mom"]),
                ("father", trio["dad"]),
            ]
            for role, pid in members:
                if not pid:
                    continue
                
                for set_id, paths in vcfs.items():
                    if pid not in set_id:
                        continue
                    for path in paths:
                        results.append({
                            "ws": label,
                            "fam": trio["fam"],
                            "sample": pid,
                            "role": role,
                            "assay": assay_type(path),
                            "type": vcf_type(path),
                            "vcf": path,
                        })
    
    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 find_vcfs_by_hpo_gregor.py HP:0001250")
        sys.exit(1)
    
    HPO = sys.argv[1]
    print(f"Looking for trios with {HPO}...", file=sys.stderr)
    
    results = find_vcfs(HPO)
    
    # header
    print("workspace\tfamily_id\tsample_id\trole\tassay\tvcf_type\tvcf_path")
    
    # sort by family then role
    role_order = {"proband": 0, "mother": 1, "father": 2}
    results.sort(key=lambda x: (x["ws"], x["fam"], role_order.get(x["role"], 9), x["assay"]))
    
    for r in results:
        print(f"{r['ws']}\t{r['fam']}\t{r['sample']}\t{r['role']}\t{r['assay']}\t{r['type']}\t{r['vcf']}")
    
    n = len(set((r["ws"], r["sample"]) for r in results))
    print(f"\nFound {len(results)} VCFs from {n} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
