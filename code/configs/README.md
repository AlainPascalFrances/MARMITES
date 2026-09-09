# Run configurations  (WP0)

One TOML file per model configuration, e.g. `lamata.toml`, `cdl.toml`. Each replaces the
command-line flags of `code/tests/run_lamata_mf6.py`; the schema and its validation live in
`code/marmites_config.py`. Machine-specific paths do **not** belong here — they go in
`code/mm_paths.py` (or the matching `MM_*` environment variables).

Section reference and the full flag-to-key map: the cookbook, §7 and Appendix A.
