# FC2026_01_02_PoP_alt_KOLF & FC2026_01_02_PoP_alt_Liver — what to commit, ignore, delete

Reference: **FC2024_11_01_PoP64_1** — everything under `02_trimming/`, `03_zUMIs/`, and the three deML files is committed; nothing in that folder is in `.gitignore`.

---

## Files actually used by `scripts/read_funnel_alt_samples.R`

| Dataset | Path / pattern | Purpose |
|--------|----------------|---------|
| **Liver** | `FC2026_01_02_PoP_alt_Liver/deML_summary_to_read_in.txt` | deML summary (used for **both** Liver and KOLF sections in the script) |
| **Liver** | `FC2026_01_02_PoP_alt_Liver/02_trimming/old.txt`, `PTO.txt` | Trimming stats |
| **Liver** | `FC2026_01_02_PoP_alt_Liver/03_zUMIs/{old,PTO}/zUMIs_output/stats/PoPalt_Liver_*.readspercell.txt` | Mapping fractions |
| **Liver** | `.../03_zUMIs/.../...kept_barcodes_binned.txt` | Barcode counts |
| **Liver** | `.../03_zUMIs/{old,PTO}/zUMIs_output/stats/PoPalt_Liver_{old,PTO}.bc.READcounts.rds` | Read counts |
| **Liver** | `.../03_zUMIs/{old,PTO}/zUMIs_output/stats/PoPalt_Liver_{old,PTO}.UMIcounts.txt` | UMI counts |
| **KOLF** | `FC2026_01_02_PoP_alt_KOLF/02_trimming/old.txt`, `PTO.txt` | Trimming stats |
| **KOLF** | `FC2026_01_02_PoP_alt_KOLF/03_zUMIs/...` | Same structure as Liver (readspercell, kept_barcodes_binned, bc.READcounts.rds, UMIcounts.txt) |

KOLF’s `deML_summary_to_read_in.txt` is **not** read by the script (both sections use the Liver deML file).

---

## What to COMMIT (match FC2024_11_01_PoP64_1)

Commit the same layout as in `data/FC2024_11_01_PoP64_1/`:

| Location | Contents |
|----------|----------|
| **Both FC2026 folders** | |
| `02_trimming/` | `old.txt`, `PTO.txt`, `cutadapt_adapt.sh` |
| `03_zUMIs/` | `BCs.txt`, `zUMIs_text.txt` |
| `03_zUMIs/old/` | `BCs_old.txt`, `*.yaml`, `*.run.yaml`, `PoPalt_*.BCstats.txt`, `*zUMIs_runlog.txt`, `*zUMIs_YAMLerror.log`, `tmp.*.SJ.out.tab` |
| `03_zUMIs/old/zUMIs_output/` | `*.BCbinning.txt`, `*kept_barcodes*.txt`, `expression/*.rds`, `expression/*.gene_names.txt`, `stats/*` (all .rds, .txt, .pdf) |
| `03_zUMIs/PTO/` | Same structure as `old/` |
| **Root of each FC2026 folder** | `deML_summary_to_read_in.txt`, `deML_summary.txt`, `deML_error.txt` |

So: everything under `02_trimming/`, `03_zUMIs/`, and the three deML files; **do not** commit `03b_genome/`.

---

## What to IGNORE (add to `.gitignore` if not already covered)

| Path / pattern | Reason |
|----------------|--------|
| `data/FC2026_01_02_PoP_alt_KOLF/03b_genome/` | STAR index (large, generated). Not matched by current `**/03b_custom_genome(s)/`. |
| `data/FC2026_01_02_PoP_alt_Liver/03b_genome/` | Same (folder doesn’t exist in Liver listing but add for consistency if you create it). |

Optional: add a single rule so any future flowcell with a STAR index under `03b_genome` is ignored, e.g. `**/03b_genome/`.

---

## What to DELETE

| Action | Item |
|--------|------|
| **Do not delete** | Any of the files listed under “COMMIT” or “Files actually used” (script would break). |
| **Optional (on disk only)** | `03b_genome/` can be removed from disk if you don’t need the index locally (rebuilt when needed); keep it ignored in git. |
| **Nothing to delete from git** | FC2026 dirs are currently untracked; there’s nothing to remove from the repo. |

---

## Summary table (commit / ignore / delete)

| Category | FC2026_01_02_PoP_alt_KOLF | FC2026_01_02_PoP_alt_Liver |
|----------|----------------------------|-----------------------------|
| **Commit** | `02_trimming/` (old.txt, PTO.txt, cutadapt_adapt.sh), `03_zUMIs/` (full tree as in FC2024_11_01), `deML_summary_to_read_in.txt`, `deML_summary.txt`, `deML_error.txt` | Same |
| **Ignore** | `03b_genome/` (STAR index) | `03b_genome/` if present |
| **Delete** | None | None |

---

## FC2024_11_01_PoP64_1 (reference)

- **Committed:** All of `02_trimming/`, `03_zUMIs/` (including yamls, logs, tmp.*.SJ.out.tab, zUMIs_output stats/expression), `deML_summary_to_read_in.txt`, `deML_summary.txt`, `deML_error.txt`. No `00_fastq`, no `03b_*`, no genome index.
- **Ignored:** Nothing under this folder; global rules apply (`**/00_fastq/*.fastq.gz`, `**/02_trimming/*.fastq.gz`, `**/03b_custom_genome(s)/`, etc.).
