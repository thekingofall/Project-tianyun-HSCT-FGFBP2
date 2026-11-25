"""Run verified saved-table plotting jobs in an isolated output directory.

All jobs read ONLY compact tables included in this GitHub repository (no multi-GB
Seurat RDS files required) and write outputs to a fresh user-specified directory.
"""
from pathlib import Path
import argparse
import json
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

EXTERNAL_JOBS = {
    "M04h": {
        "module": "GSE216187_effector_DE",
        "script": "14_heatmap_markers16_compact_v4.R",
        "description": "Fig. 4H | GSE216187 16-marker CR vs REL CD8 effector heatmap (Z-scores)",
        "mode": "archived",
    },
    "M04i": {
        "module": "GSE222633_CR_REL",
        "script": "02_redraw_CR_REL.R",
        "description": "Fig. 4I | GSE222633 patient-level longitudinal FGFBP2/ADGRG1 trajectories (CR vs REL)",
        "mode": "archived",
    },
    "S05b": {
        "module": "GSE216187_markers",
        "script": "11_render_marker_heatmap_final.R",
        "description": "Fig. S5B | GSE216187 24-celltype canonical marker expression heatmap",
        "mode": "archived",
    },
    "S05c": {
        "module": "GSE216187_similarity",
        "script": "33_direct_full15_plot.R",
        "description": "Fig. S5C | GSE216187 Spearman transcriptional similarity heatmap vs primary CD8 subsets",
        "mode": "wrap_ggplot",
        "out_name": "S05c_GSE216187_full15_similarity_heatmap",
        "width": 9.2,
        "height": 7.4,
    },
    "S05ef": {
        "module": "GSE222633_CD8",
        "script": "21_plot_FCD8.R",
        "description": "Fig. S5E-F | GSE222633 CD8 UMAP and canonical marker dot plot",
        "mode": "archived",
    },
    "S05f": {
        "module": "GSE222633_CD8",
        "script": "25_marker_dotplot_v2.R",
        "description": "Fig. S5F | GSE222633 CD8 revised canonical marker dot plot",
        "mode": "archived",
    },
    "S05g": {
        "module": "GSE222633_longitudinal",
        "script": "06_redraw_all_figures.R",
        "description": "Fig. S5G | GSE222633 longitudinal Teff FGFBP2/ADGRG1 expression and functional scores",
        "mode": "archived",
    },
}

PRIMARY_JOBS = {
    "M01": "Fig. 1B-E | Primary PBMC UMAP, canonical marker dot plot, lineage composition, and stage density",
    "M02": "Fig. 2A-D & S2D | Primary T cell 15-subset UMAP, marker dot plot, composition, and functional scores",
    "M03": "Fig. 3A-E,G-K & S3A-D | Primary CD8/CD4 trajectory embeddings and state compositions",
    "S06e": "Fig. S6E | Inhibitory receptor and TOX positivity percentages across T cell states",
    "primary_all": "Run all primary saved-table verification modules (M01, M02, M03, S06e)",
}


def sanitize_archived_r(code_text: str, target_dir: Path, spec: dict) -> str:
    """Make archived server R scripts portable inside an isolated output folder."""
    s = code_text
    # Remove server-specific fontconfig and setwd calls
    s = re.sub(r"^.*Sys\.setenv\(FONTCONFIG_FILE=.*\n", "", s, flags=re.M)
    s = re.sub(r"^.*read_json\('resources/figure_runtime\.json'.*\n", "", s, flags=re.M)
    s = re.sub(r"^setwd\([^\n]+\)\s*$", "", s, flags=re.M)
    s = re.sub(r"^path\s*<-\s*['\"][^'\"]+['\"]", "path <- getwd()", s, flags=re.M)
    # Prevent non-interactive Rscript from opening default PostScript Rplots.pdf on bare print(p)
    s = re.sub(r"^\s*print\(p\)\s*$", "", s, flags=re.M)
    # Ensure ggsave PNG calls explicitly pass device=grDevices::png (without touching grDevices::png calls)
    fixed_lines = []
    for line in s.splitlines():
        if "ggsave(" in line and "grDevices::png(" not in line:
            line = re.sub(
                r"(?<!device=grDevices::png,)type\s*=\s*(['\"])cairo\1",
                r"device=grDevices::png,type=\1cairo\1",
                line,
            )
        fixed_lines.append(line)
    s = "\n".join(fixed_lines) + ("\n" if s.endswith("\n") else "")
    if spec["mode"] == "wrap_ggplot":
        # Strip trailing bare `p` print statement that triggers Rplots.pdf
        s = re.sub(r"^\s*p\s*$", "", s, flags=re.M)
        out_name = spec["out_name"]
        w, h = spec["width"], spec["height"]
        s += (
            f"\npng_type <- if (isTRUE(capabilities('cairo'))) 'cairo' else 'Xlib'\n"
            f"grDevices::png('figures/png/{out_name}.png', width={w}, height={h}, units='in', res=300, type=png_type, bg='white'); print(p); grDevices::dev.off()\n"
            f"grDevices::cairo_pdf('figures/pdf/{out_name}.pdf', width={w}, height={h}, bg='white'); print(p); grDevices::dev.off()\n"
            f"cat('WRAPPED_GGPLOT_READY: {out_name}\\n')\n"
        )
    # Ensure resources/figure_runtime.json exists if any script still checks it
    res_dir = target_dir / "resources"
    if res_dir.is_dir() and not (res_dir / "figure_runtime.json").exists():
        (res_dir / "figure_runtime.json").write_text(json.dumps({"fontconfig": ""}), encoding="utf-8")
    return s


def run_job(panel: str, output: Path, rscript: str, prepare_only: bool) -> None:
    target = output.resolve()
    if target.exists():
        raise SystemExit(f"Output directory already exists: {target}. Please choose a new directory.")

    if panel in PRIMARY_JOBS:
        target.mkdir(parents=True, exist_ok=True)
        mod_arg = "all" if panel == "primary_all" else panel
        script = ROOT / "00_Settings" / "03_reproduce_primary_panels.R"
        cmd = [rscript, str(script), mod_arg, str(ROOT), str(target)]
        if prepare_only:
            print("Prepared primary command:", " ".join(cmd))
            return
        subprocess.run(cmd, cwd=target, check=True)
        print("Results:", target / "figures")
        return

    spec = EXTERNAL_JOBS[panel]
    inputs = ROOT / "00_Global_Data" / "Shared_Inputs" / spec["module"]
    code = ROOT / "Main_Figure_04" / "00_Preprocess" / spec["module"] / spec["script"]
    shutil.copytree(inputs, target)
    for d in ["figures/png", "figures/pdf", "logs"]:
        (target / d).mkdir(parents=True, exist_ok=True)
    runfile = target / "run_plot.R"
    runfile.write_text(sanitize_archived_r(code.read_text(encoding="utf-8-sig"), target, spec), encoding="utf-8")
    if prepare_only:
        print("Prepared:", runfile)
        return
    subprocess.run([rscript, str(runfile)], cwd=target, check=True)
    print("Results:", target / "figures")


def main() -> None:
    all_choices = list(PRIMARY_JOBS.keys()) + list(EXTERNAL_JOBS.keys())
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", nargs="?", choices=all_choices, help="Panel or module ID to reproduce")
    p.add_argument("--list", action="store_true", help="List all supported saved-table reproduction jobs")
    p.add_argument("--output", help="Path to a new output directory")
    p.add_argument("--rscript", default="Rscript", help="Path to Rscript executable (default: Rscript)")
    p.add_argument("--prepare-only", action="store_true", help="Stage inputs and portable script without executing R")
    a = p.parse_args()

    if a.list or not a.panel:
        print("Primary cohort saved-table verification jobs:")
        for k, desc in PRIMARY_JOBS.items():
            print(f"  {k:<12} {desc}")
        print("\nExternal validation cohort saved-table jobs:")
        for k, spec in EXTERNAL_JOBS.items():
            print(f"  {k:<12} {spec['description']}")
        return

    if not a.output:
        p.error("--output is required when running a panel job.")
    run_job(a.panel, Path(a.output), a.rscript, a.prepare_only)


if __name__ == "__main__":
    main()
