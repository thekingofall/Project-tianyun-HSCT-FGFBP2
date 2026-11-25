"""Run selected archived table-based plots in a new output directory."""
from pathlib import Path
import argparse, shutil, subprocess, re
ROOT=Path(__file__).resolve().parents[1]
JOBS={
 'M04h':('GSE216187_effector_DE','14_heatmap_markers16_compact_v4.R'),
 'S05ef':('GSE222633_CD8','21_plot_FCD8.R'),
 'S05f':('GSE222633_CD8','25_marker_dotplot_v2.R'),
}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('panel',choices=JOBS);p.add_argument('--output',required=True);p.add_argument('--rscript',default='Rscript');p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
 module,script=JOBS[a.panel];target=Path(a.output).resolve()
 if target.exists():p.error('Output directory already exists; choose a new directory.')
 inputs=ROOT/'00_Global_Data/Shared_Inputs'/module
 code=ROOT/'Main_Figure_04/00_Preprocess'/module/script
 shutil.copytree(inputs,target)
 s=code.read_text(encoding='utf-8-sig')
 s=re.sub(r'^.*Sys\.setenv\(FONTCONFIG_FILE=.*\n','',s,flags=re.M)
 s=re.sub(r'^setwd\([^\n]+\)\s*$', '',s,flags=re.M)
 s=re.sub(r"^path\s*<-\s*['\"][^'\"]+['\"]",'path <- getwd()',s,flags=re.M)
 for d in ['figures/png','figures/pdf','logs']:(target/d).mkdir(parents=True,exist_ok=True)
 runfile=target/'run_plot.R';runfile.write_text(s,encoding='utf-8')
 if a.prepare_only:print('Prepared:',runfile);return
 subprocess.run([a.rscript,str(runfile)],cwd=target,check=True)
 print('Results:',target/'figures')
if __name__=='__main__':main()
