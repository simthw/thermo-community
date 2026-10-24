# Simulation Registry

`runs.toml` is the source of truth for named groups of simulation outputs.

Use it for any figure script that compares multiple runs. The raw data can
remain in the `simulations/` tree, but the question "which dirs
belong to this analysis?" should be answered here instead of in comments or
hard-coded arrays inside plotting scripts.

## Directory Policy

- Raw simulation outputs: timestamped run directories with `config.toml`,
  `summary.csv`, `paras/`, and optionally `trajectories/`.
- Derived analysis/figure outputs: paths listed under each run set's `outputs`
  table, below `results/<run_set_name>/`.
- Presentation exports: keep under `results/presentation/`.
- Older runs left in place: do not register them unless they are current-model
  output; the registry holds the datasets the project actually uses.

## Adding A New Run Set

Add one `[sets.<name>]` table with:

- `description`: one sentence stating what the data are.
- `kind`: for example `simulation_sweep`, `ablation_sweep`, or
  `two_condition_comparison`.
- `status`: `candidate` or `canonical`.
- `thermodynamic_feedback`: explicit `true` or `false`; never infer the mode from a missing field.
- `paired_seeds = true` when every condition must contain the same seed set;
  before promotion, also run `validate_paired_run_directories` to verify that
  each shared seed has identical reaction and strain parameters.
- `varying_config_keys`: config fields intentionally allowed to differ across
  paired conditions; use `["sulfate_mg_L"]` for the standard sulfate sweep.
- `dirs`: ordered list of `{ label, so4_mg_L, path }` entries.
- optional `[sets.<name>.outputs]`: derived CSV/figure directories for scripts.

Prefer names such as `full_ensemble_2026_08_01` over bare timestamps: the name
says what the data mean, and the raw directory path records when it was made.
New sets start as `candidate` and become `canonical` only after the gates in
`../../EXECUTION_RUNBOOK.md` pass.
