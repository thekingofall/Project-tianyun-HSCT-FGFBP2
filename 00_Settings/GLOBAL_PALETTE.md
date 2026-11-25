# Global Palette Registry / 全局配色锁定表

All label-to-HEX mappings below are locked directly to the submitted manuscript figures and archived plotting scripts. No categorical color was reassigned during repository organization.

- **Machine-readable table**: [global_palette.tsv](global_palette.tsv)
- **Cohort-specific registries**: preserved verbatim under `00_Global_Data/Shared_Inputs/*/palette_registry.tsv`.

---

## 1. Primary Cohort: PBMC Lineages (`Fig. 1B, 1D, 1E; Fig. S1B`)

| Order | Lineage | HEX |
|---:|---|---|
| 1 | `HSC/MPPs` | `#49B2BE` |
| 2 | `T` | `#699ECA` |
| 3 | `NK` | `#0098E2` |
| 4 | `B` | `#FB5039` |
| 5 | `Monocytes` | `#C07DD9` |
| 6 | `Granulocytes` | `#D97DAA` |
| 7 | `cDC` | `#FF8C00` |
| 8 | `pDC` | `#DE3163` |
| 9 | `Platelets` | `#8582BD` |

## 2. Primary Cohort: T Cell Subsets (`Fig. 2A, 2C; Fig. 3A, 3D, 3G, 3J; Fig. S2B`)

| Order | T Cell Subset | HEX |
|---:|---|---|
| 1 | `CD4 Naive` | `#A3CFE8` |
| 2 | `CD4 Tscm` | `#8AC4C1` |
| 3 | `CD4 Tcm` | `#41A6AB` |
| 4 | `CD4 Tem` | `#DFAA81` |
| 5 | `CD4 Teff` | `#CA8652` |
| 6 | `CD4 Treg` | `#246B9A` |
| 7 | `CD8 Naive` | `#CDB4D8` |
| 8 | `CD8 Tcm` | `#EC8EC8` |
| 9 | `CD8 Tem` | `#764B88` |
| 10 | `CD8 Teff` | `#A0478D` |
| 11 | `CD8 Tex` | `#F7B7A3` |
| 12 | `Trm-like` | `#D48B8B` |
| 13 | `NKT` | `#F37A6A` |
| 14 | `MAIT` | `#9DAE55` |
| 15 | `γδT` | `#64AE55` |

## 3. Primary Cohort: Sampling Stages & Trajectory States (`Fig. 2D–F; Fig. 3C, 3I; Fig. S1–S3`)

| Group | Label | HEX |
|---|---|---|
| Sampling Stage | `Pre-HSCT` (`Pre`) | `#7030A0` |
| Sampling Stage | `HD` | `#0070C0` |
| Sampling Stage | `Post-1M` | `#FF5050` |
| Sampling Stage | `Post-3M` | `#F78832` |
| Sampling Stage | `Post-6M` | `#70AD47` |
| Trajectory State | `S1` | `#93A9FB` |
| Trajectory State | `S2` | `#E285E6` |
| Trajectory State | `S3` | `#6FBA8A` |

## 4. External Validation Cohort: GSE216187 (`Fig. 4H; Fig. S5A–D`)

| Group | Label | HEX |
|---|---|---|
| CD8 Compartment | `CD8 NV` / `CD8 hobit` / `CD8 mem. 1–3` / `CD8 eff. 1–2` / `CD8 IFN` | `#DFC4DE`, `#B092B2`, `#D475A8`, `#A8C5DB`, `#7EB5D6`, `#569ECF`, `#3D7AB8`, `#6B429E` |
| CD4 & Unconventional T | `CD4 NV` / `CD4 mem.` / `Th17` / `Treg` / `CD4 IFN` / `γδT` / `MAIT` | `#E38888`, `#E36442`, `#D9382B`, `#A81C1C`, `#7A1122`, `#8B52A3`, `#5C2D82` |
| HSPC & Myeloid | `MLP` / `MEP` / `NP` / `MDP` / `pDCs` / `cDCs` / `preB` / `pro/pre T` / `MP/mono` | `#1F8A8A`, `#229EAB`, `#3BC4BE`, `#42B883`, `#7BC993`, `#4DB86E`, `#2E9E52`, `#1B7D3A`, `#115C27` |
| Outcome (`Fig. 4H`) | `CR` / `REL` | `#B15A50` / `#61749B` |
| Z-score Scale (`Fig. 4H`) | `-2` / `0` / `+2` | `#5C53B6` / `#F7F7F7` / `#F16643` |

## 5. External Validation Cohort: GSE222633 (`Fig. 4I; Fig. S5E–G`)

| Group | Label | HEX |
|---|---|---|
| CD8 Subsets (`Fig. S5E–F`) | `CD8 Naive` / `CD8 Tem` / `CD8 Teff` / `CD8 Cycling` / `CD8 Activated` | `#DFC4DE`, `#7EB5D6`, `#3D7AB8`, `#C69A45`, `#A34F65` |
| Unconventional / Review | `MAIT-like T` / `gdT-like T` / `CD8 T (low RNA)` / `CD4-like T (review)` / `NK/T-like (review)` | `#5C2D82`, `#8B52A3`, `#BCC6CC`, `#E38888`, `#5C6D73` |
| Donor–Recipient Pairs | `Pair A` / `Pair B` / `Pair C` / `Pair D (REL)` / `Pair E` | `#3B6C80`, `#B47A4E`, `#68875A`, `#A34F65`, `#7C6A98` |
| Clinical Outcome | `CR` / `REL` | `#2878B5` / `#D64B40` |
