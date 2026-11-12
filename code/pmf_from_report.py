"""Parse pasted report_pmf.py output (all four systems) into pmf.json.
Usage: python3 pmf_from_report.py report.txt pmf.json"""
import sys,re,json
txt=open(sys.argv[1],encoding="utf-8").read(); P={}; cur=None
for ln in txt.splitlines():
    m=re.match(r"\s*(PTFE|PVDF|ECTFE|PE)(?:_PFOA)?(?:_\d+_\d+)?\s*$",ln)
    if m: cur=m.group(1); P[cur]={"ns":5}; continue
    if cur is None: continue
    m=re.search(r"dG_ads\s*=\s*(-?\d+\.\d+)\s*\+/-\s*(\d+\.\d+)",ln)
    if m: P[cur]["dG"]=float(m.group(1)); P[cur]["err"]=float(m.group(2))
    m=re.search(r"인접 창 중첩.*최소\s*(\d+)%",ln)
    if m: P[cur]["overlap"]=float(m.group(1))
    m=re.search(r"(\d+)\s*ns",ln)
    if m and "연장" not in ln and "sampl" in ln.lower(): P[cur]["ns"]=int(m.group(1))
missing=[p for p in ("PTFE","PVDF","ECTFE","PE") if p not in P or "dG" not in P[p]]
if missing: print("WARNING missing:",missing)
json.dump(P,open(sys.argv[2],"w"),indent=1); print(json.dumps(P,indent=1))
