import subprocess
import sys
import time
from pathlib import Path

from _bootstrap import pd, setup, write_table


def main():
    args, _, out = setup("Run all source-data analyses", "99_run_log")
    script_root = Path(__file__).resolve().parent
    scripts = [script_root / "00_validate_inputs.py"]
    scripts += [script_root / f"{number:02d}_{name}" for number, name in [
        (1, "cohort_clinical_phenotypes.py"),
        (2, "major_cell_atlas_composition.py"),
        (3, "cime_programs.py"),
        (4, "bcell_states_differentiation.py"),
        (5, "bcr_repertoire.py"),
        (6, "cd4_hiv_transcript_host_response.py"),
        (7, "cd4_tcr_repertoire.py"),
        (8, "cd4_cross_compartment.py"),
        (9, "cd8_state_composition.py"),
        (10, "cd8_tcr_overview_sankey.py"),
        (11, "cd8_tcr_recurrence_similarity.py"),
        (12, "cellchat_global_networks.py"),
        (13, "myeloid_microglia_states.py"),
        (14, "microglia_cellchat_differential.py"),
        (15, "proteomics_differential_shared_core.py"),
        (16, "proteomics_nmf_modules.py"),
        (17, "bcr_proteome_integration.py"),
        (18, "proteomic_diagnostic_models.py"),
        (19, "cime_proteome_symptom_triad.py"),
        (20, "statistical_reanalysis.py"),
    ]]
    rows = []
    panels = []
    for script in scripts:
        started = time.time()
        command = [sys.executable, str(script), "--data-root", str(args.data_root), "--output-root", str(args.output_root), "--seed", str(args.seed)]
        completed = subprocess.run(command, capture_output=True, text=True)
        rows.append({"script": script.name, "return_code": completed.returncode, "seconds": time.time() - started, "stderr": completed.stderr.strip()})
        write_table(pd.DataFrame(rows), out / "run_log.tsv")
        if completed.returncode:
            raise RuntimeError(f"{script.name} failed; see run_log.tsv")
        for figure in sorted(args.output_root.rglob("*.pdf")):
            relative = figure.relative_to(args.output_root)
            if not any(row["figure"] == str(relative) for row in panels):
                panels.append({"script": script.name, "figure": str(relative)})
    write_table(pd.DataFrame(panels), out / "panel_manifest.tsv")


if __name__ == "__main__":
    main()
