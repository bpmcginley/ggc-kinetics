# Data attribution

This package redistributes **no third-party data files**. The numbers below are typed into the tests or the
benchmark, with their sources.

| Used in | Source | Licence | What is used |
|---|---|---|---|
| tests (invertase) | Michaelis & Menten 1913, Table 1 as translated in Johnson & Goody 2011, Biochemistry, PMC3381512 (doi:10.1021/bi201284u) | public domain (1913 data) | 31 time / fraction points |
| tests (PON1) | Goličnik & Bavec 2020, J. Enzyme Inhib. Med. Chem., PMC6896510 (doi:10.1080/14756366.2019.1695792), Table 1 | CC BY 4.0 | 12 published % error values |
| benchmark only (not shipped) | ICEKAT supplements (Olp, Kalous & Smith 2020, BMC Bioinformatics, PMC7222511, doi:10.1186/s12859-020-3513-y) and the ICEKAT `test.csv` example | CC BY 4.0 | parsed and fitted locally; see `benchmarks/REAL_DATA_BENCHMARKS.md` |
| benchmark only (not shipped) | HRP / ABTS Source Data of Kont et al. 2020, Nat. Commun. (PMC7666214, doi:10.1038/s41467-020-19561-8) | CC BY 4.0 | plain MM fits per series, locally |
| benchmark only (not shipped) | DaRUS datasets doi:10.18419/darus-3867 (Häussler 2024, SLAC / ABTS) and doi:10.18419/darus-2467 (Ngubane 2022, laccase / ABTS) | CC BY 4.0 | laccase / ABTS plates fitted locally |

`tests/data/` holds MODELED (simulated) plate files from GGC-inf Cycle 33 and the reports Cycle 34 stored for them.
They are not measurements.

If you redistribute any CC BY data next to this package, cite the original authors and DOIs above and state the
licence.
