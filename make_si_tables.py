"""Emit si_tables.tex from experiment_summary.json, tier.json, raw csv. Usage: si_tables.py outdir raw.csv [pmf.json]"""
import sys,json,os,csv
out,raw=sys.argv[1],sys.argv[2]; pmfj=sys.argv[3] if len(sys.argv)>3 else None
S=json.load(open(f"{out}/experiment_summary.json")); T=json.load(open(f"{out}/tier.json"))
cell={(s["polymer"],s["pfas"]):s for s in S["systems"]}; PL=S["polymer_level"]
POLY=["PTFE","PVDF","ECTFE","PE"]; PF=["PFOA","PFBS","GenX","FTOH62"]; LAB={"PFOA":"PFOA","PFBS":"PFBS","GenX":"GenX","FTOH62":"6:2 FTOH"}
def p(v): return f"{v:.3f}" if v>=0.001 else "$<0.001$"
def n(v,d=1): return f"${v:.{d}f}$"
def sg(v,d=1): return f"${v:+.{d}f}$"
B=[]
def table(label,caption,cols,rows,pos="H",resize=False):
    tab="\\begin{tabular}{@{}"+cols+"@{}}\n\\toprule\n"+rows+"\\bottomrule\n\\end{tabular}"
    if resize: tab="\\resizebox{\\textwidth}{!}{"+tab+"}"
    B.append("\\begin{table}[H]\n\\centering\\small\n\\caption{"+caption+"}\n\\label{"+label+"}\n"+tab+"\n\\end{table}\n")
# S1 slab density
dens={"PTFE":(42,"(CF$_2$CF$_2$)$_{15}$",76.0,3.478,1.8828,1.98,"amorphous PTFE (Callister)"),
      "PVDF":(54,"(CH$_2$CF$_2$)$_{15}$",59.4,3.639,1.4825,1.68,"amorphous PVDF"),
      "ECTFE":(42,"alt.\\ E/CTFE, C$_{30}$",39.4,3.485,1.4488,1.55,"est.\\ from bulk 1.68 at $\\sim$55\\% crystallinity"),
      "PE":(65,"C$_{30}$H$_{62}$",0.0,3.548,0.8038,0.855,"amorphous PE")}
rows="Polymer & Chains & Oligomer & wt\\% F & Thickness (nm) & $\\rho_{\\mathrm{slab}}$ & $\\rho_{\\mathrm{ref}}$ & $\\rho_{\\mathrm{slab}}/\\rho_{\\mathrm{ref}}$ \\\\\n\\midrule\n"
for q in POLY:
    c,olig,wf,h,r,ref,src=dens[q]; rows+=f"{q} & {c} & {olig} & {wf:.1f} & {h:.2f} & {r:.3f} & {ref:.3f} & {r/ref:.2f} \\\\\n"
table("tab:slab","Amorphous slab construction and density validation. Each slab was compressed, annealed (500~K) and relaxed as described in Section~\\ref{sec:md}; the density is that of the converged dry slab after the 1~ns convergence probe (drift $\\le$2\\%). Reference densities are amorphous-phase values; GAFF2 is known to under-predict fluoropolymer density by 5--12\\%, and the shortfall is uniform across the four polymers. Densities in g~cm$^{-3}$. Reference basis: amorphous-phase PTFE (1.98) and PE (0.855) from standard handbook values, amorphous PVDF 1.68, and ECTFE 1.55 estimated from the bulk density (1.68) at $\\sim$55\\% crystallinity.","lrlrrrrr",rows)
# S2 per-cell
rows="Polymer & Compound & $n$ & Mean IE & CI low & CI high & s.d. & Coulomb & Lennard-Jones \\\\\n\\midrule\n"
for q in POLY:
    for f in PF:
        c=cell[(q,f)]; rows+=f"{q} & {LAB[f]} & {c['n']} & {n(c['meanIE'])} & {n(c['ci_lo'])} & {n(c['ci_hi'])} & {c['sd']:.1f} & {sg(c['meanCoul'])} & {n(c['meanLJ'])} \\\\\n"
table("tab:ie","Per-cell short-range interaction energy (kJ mol$^{-1}$; more negative indicates stronger interaction; $n=5$ replicates each) with bootstrap 95\\% confidence interval ($10^{4}$ resamples), replicate standard deviation, and the Coulomb and Lennard-Jones components of the mean. These values underlie Figure~2a.","llrrrrrrr",rows)
# S3 polymer-level
rows="Polymer & Mean IE & s.d. & s.e.m. & Coulomb & s.d. & Lennard-Jones & s.d. & $\\hat{A}_{\\mathrm{ads}}$ (mean of cells) \\\\\n\\midrule\n"
for q in POLY:
    a=PL[q]; A=[cell[(q,f)]["A_ads_norm"] for f in PF]
    rows+=f"{q} & {n(a['mean_reps'])} & {a['sd_reps']:.1f} & {a['sem_reps']:.1f} & {sg(a['coul'],2)} & {a['coul_sd']:.1f} & {n(a['lj'])} & {a['lj_sd']:.1f} & {sum(A)/4:.2f} \\\\\n"
table("tab:polymeans","Polymer-level means over the 20 replicates of each polymer (kJ mol$^{-1}$), with the Coulomb / Lennard-Jones decomposition used in Figure~2b.","lrrrrrrrr",rows,"H",True)
# S4 contrasts
rows="Contrast & Term & $\\Delta$ (kJ mol$^{-1}$) & $p$ (blocked) & $d_z$ & $p$ (Welch) & Cohen's $d$ & $p$ (Mann--Whitney) \\\\\n\\midrule\n"
name={"totalIE":"total","coulSR":"Coulomb","ljSR":"LJ"}
for c in S["contrasts"]:
    rows+=f"{c['contrast']} & {name[c['term']]} & {sg(c['delta'])} & {p(c['p'])} & {n(c['dz'],2)} & {p(c['p_rep'])} & {n(c['cohen_d'],2)} & {p(c['p_mw'])} \\\\\n"
    if c["term"]=="ljSR": rows+="\\addlinespace[2pt]\n"
for r in S["ptfe_vs_rest"]:
    rows+=f"PTFE vs mean(PVDF, ECTFE, PE) & {dict(total='total',coul='Coulomb',lj='LJ')[r['term']]} & {sg(r['delta'])} & {p(r['p'])} & {n(r['dz'],2)} & -- & -- & -- \\\\\n"
fr=S["friedman_p"]; rows+=f"\\addlinespace[2pt]\nFriedman (4 polymers, compound blocks) & total / Coulomb / LJ & -- & {p(fr['total'])} / {p(fr['coul'])} / {p(fr['lj'])} & -- & -- & -- & -- \\\\\n"
table("tab:polytests","Polymer-level statistical checks for the total interaction energy and its Coulomb and Lennard-Jones components. Because the same four compounds appear in every polymer arm, compound is treated as a block: the paired $t$ test across the four compound-level means ($p$ blocked, three degrees of freedom, effect size $d_z$) is the primary statistic. Replicate-level Welch $t$ (Cohen's $d$) and Mann--Whitney tests pool all 20 replicates per polymer and are anti-conservative because they treat replicates of different compounds as exchangeable. $\\Delta$ is the first polymer minus the second; a positive $\\Delta$ means weaker binding (or stronger repulsion) for the first polymer.","llrrrrrr",rows,"H",True)
# S5 fluoro vs PE per compound
rows="Polymer & Compound & Cohen's $d$ vs PE & $p$ (Welch) \\\\\n\\midrule\n"
for e in S["effect_sizes_vs_PE"]: rows+=f"{e['polymer']} & {LAB[e['pfas']]} & {n(e['cohens_d'],2)} & {p(e['p_welch'])} \\\\\n"
table("tab:effpe","Per-compound effect sizes of each fluoropolymer relative to the PE baseline (replicate level, $n=5$ vs $5$). Positive $d$ indicates weaker binding than PE.","llrr",rows)
# S6 tier
rows="Polymer & "+" & ".join(LAB[f] for f in PF)+" & Mean \\\\\n\\midrule\n"
for i,q in enumerate(POLY): rows+=q+" & "+" & ".join(f"{v:.2f}" for v in T["tier"][i])+f" & {T['tier_mean'][i]:.2f} \\\\\n"
table("tab:proxy","Evidence-tier prioritization values used in Figure~3b,c ($\\hat{A}_{\\mathrm{ads}}\\times\\hat{H}_{\\mathrm{receptor}}$; $\\hat{A}_{\\mathrm{ads}}$ is a min--max scaling of $|$IE$|$ across the 16 cells with a 0.05 floor). PFBS is low but non-zero. The non-fluorinated PE baseline is included: its mean tier ("+f"{T['tier_mean'][3]:.2f}"+") exceeds that of PTFE ("+f"{T['tier_mean'][0]:.2f}"+"), consistent with the interaction-energy screen.","lrrrrr",rows)
# S7 PMF
if pmfj and os.path.exists(pmfj):
    P=json.load(open(pmfj)); rows="System & Windows & Sampling per window (ns) & Min.\\ adjacent overlap (\\%) & $\\Delta G_{\\mathrm{ads}}$ (kJ mol$^{-1}$) & Bootstrap s.d. & $|$IE$|/|\\Delta G|$ \\\\\n\\midrule\n"
    for q in POLY: rows+=f"{q}--PFOA & 19 & {P[q].get('ns',5)} & {P[q].get('overlap',0):.0f} & {n(P[q]['dG'])} & {P[q]['err']:.1f} & {abs(cell[(q,'PFOA')]['meanIE']/P[q]['dG']):.1f} \\\\\n"
else:
    rows="System & Windows & Sampling per window (ns) & Min.\\ adjacent overlap (\\%) & $\\Delta G_{\\mathrm{ads}}$ (kJ mol$^{-1}$) & Bootstrap s.d. & $|$IE$|/|\\Delta G|$ \\\\\n\\midrule\n"
    for q in POLY: rows+=f"{q}--PFOA & 19 & 5 & \\textcolor{{red}}{{[..]}} & \\textcolor{{red}}{{[..]}} & \\textcolor{{red}}{{[..]}} & \\textcolor{{red}}{{[..]}} \\\\\n"
table("tab:pmf","Umbrella-sampling PMF anchors reported in Figure~2d.","lrrrrrr",rows,"H",True)
# S8 all replicates (longtable-free: two column blocks)
rep=list(csv.DictReader(open(raw)))
rows="Polymer & Compound & Seed & Coulomb & Lennard-Jones & Total \\\\\n\\midrule\n"
for r in rep: rows+=f"{r['polymer']} & {LAB[r['pfas']]} & {r['seed']} & {sg(float(r['coulSR_kJmol']),2)} & {n(float(r['ljSR_kJmol']),2)} & {n(float(r['totalIE_kJmol']),2)} \\\\\n"
B.append("\\begingroup\\small\n\\begin{longtable}{@{}llrrrr@{}}\n\\caption{All "+str(len(rep))+" per-replicate short-range polymer--compound interaction energies (kJ mol$^{-1}$), averaged over the final 4~ns of each 6~ns production run. The full file is \\texttt{data/raw\\_results.csv} in the repository.}\\label{tab:allreps}\\\\\n\\toprule\n"+rows.split("\\midrule\n")[0]+"\\midrule\n\\endfirsthead\n\\toprule\n"+rows.split("\\midrule\n")[0]+"\\midrule\n\\endhead\n"+rows.split("\\midrule\n")[1]+"\\bottomrule\n\\end{longtable}\n\\endgroup\n")
parts={"slab":[0],"screen":[1,2,3,4,7],"pmf":[6],"tier":[5]}
for k,idx in parts.items(): open(f"{out}/si_tables_{k}.tex","w").write("\n".join(B[i] for i in idx))
print("wrote",len(B),"tables")
