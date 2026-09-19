"""
export_latex_table.py
Formats the full 56-system evaluation results into a LaTeX table for main.tex
"""
import pandas as pd

df = pd.read_csv("data/kepler_population_results.csv")

with open("data/latex_catalog_table.tex", "w", encoding="utf-8") as f:
    f.write("\\begin{longtable}{r l l r r r r r r l l}\n")
    f.write("\\caption{Complete population evaluation catalog ($N=56$) spanning four decades of transit depth. $\\text{SNR}_{\\text{trad}}$ and $\\text{VESTA}$ denote traditional matched-filter SNR and variance-stabilized score, respectively. $\\mathcal{P}_{\\text{corr}}$ is the shape-corrected logarithmic penalty, and Tier indicates the calibrated statistical category based on the null $F$-distribution.}\\label{tab:full_population_catalog}\\\\\n")
    f.write("\\hline\\hline\n")
    f.write("\\textbf{\\#} & \\textbf{Target} & \\textbf{Host Star} & \\textbf{Depth (ppm)} & \\textbf{$P$ (d)} & \\textbf{$T_{\\text{dur}}$ (h)} & \\textbf{T-SNR} & \\textbf{VESTA} & \\textbf{$\\mathcal{P}_{\\text{corr}}$} & \\textbf{Tier} & \\textbf{Class} \\\\\n")
    f.write("\\hline\n")
    f.write("\\endfirsthead\n")
    f.write("\\hline\\hline\n")
    f.write("\\textbf{\\#} & \\textbf{Target} & \\textbf{Host Star} & \\textbf{Depth (ppm)} & \\textbf{$P$ (d)} & \\textbf{$T_{\\text{dur}}$ (h)} & \\textbf{T-SNR} & \\textbf{VESTA} & \\textbf{$\\mathcal{P}_{\\text{corr}}$} & \\textbf{Tier} & \\textbf{Class} \\\\\n")
    f.write("\\hline\n")
    f.write("\\endhead\n")
    f.write("\\hline\n")
    f.write("\\endfoot\n")

    current_stratum = ""
    for _, r in df.iterrows():
        strat = r['depth_stratum']
        if strat != current_stratum:
            current_stratum = strat
            f.write(f"\\hline\n\\multicolumn{{11}}{{c}}{{\\textit{{{strat}}}}} \\\\\n\\hline\n")
        
        rank = int(r['rank'])
        target = r['target_name']
        star = r['star_name']
        depth = f"{r['literature_depth_ppm']:,.1f}"
        period = f"{r['period']:.3f}"
        dur = f"{r['duration_hours']:.1f}"
        tsnr = f"{r['t_snr']:.1f}"
        vscore = f"{r['log_snras']:.1f}"
        pen = f"{r['penalty_corr']:.4f}"
        tier = r['tier'].replace(" (Clean)", "").replace(" (Review)", "").replace(" (Veto)", "")
        tclass = "Planet" if r['label'] == 1 else "Artifact"

        f.write(f"{rank:2d} & {target:<16} & {star:<12} & {depth:>9} & {period:>7} & {dur:>5} & {tsnr:>6} & {vscore:>6} & {pen:>7} & {tier:<7} & {tclass} \\\\\n")

    f.write("\\hline\\hline\n")
    f.write("\\end{longtable}\n")

print("Generated data/latex_catalog_table.tex successfully!")
