# Shared run-set registry helpers.

# Include from scripts:
#   include(joinpath(@__DIR__, "..", "simulation_registry.jl"))

using CSV
using DataFrames
using Dates
using JLD2
using Logging
using SHA
using TOML

# Parameter JLD2 files embed the full-simulation module types (e.g.
# Analysis.AssemblyCook.FullParameters), which are not loaded in light figure
# scripts. JLD2 then warns once per type about reconstructing an equivalent
# anonymous struct. The reconstruction preserves all fields, so the loads are
# semantically identical; we silence the cosmetic warning here.
function load_parameters_silently(path::AbstractString)
  return with_logger(NullLogger()) do
    JLD2.load(path, "ps")
  end
end

const SIMULATION_PROJECT_ROOT = rstrip(normpath(joinpath(@__DIR__, "..")), ['\\', '/'])
const SIMULATION_REGISTRY_PATH = joinpath(SIMULATION_PROJECT_ROOT, "simulations", "0_rule", "runs.toml")

function file_sha256(path::AbstractString)
  isfile(path) || error("Cannot hash missing file: $path")
  return bytes2hex(SHA.sha256(read(path)))
end

function manual_source_id(run_dirs::Vector{String})
  normalized = abspath.(normpath.(run_dirs))
  digest = bytes2hex(SHA.sha256(Vector{UInt8}(codeunits(join(normalized, "\n")))))[1:10]
  return join(basename.(normalized), "__") * "__" * digest
end

"""
    infer_run_set_from_dirs(dirs) -> Union{String,Nothing}

Find the registered run set whose directory list contains ALL the given
directories. Used by figure scripts so that manual directory arguments still
resolve to a named run set (and thus output under `results/<run_set>/`)
instead of a hash-based manual source id. Returns `nothing` if no single run
set contains all directories.
"""
function infer_run_set_from_dirs(dirs::Vector{String})
  abs_dirs = abspath.(normpath.(dirs))
  for name in list_run_sets()
    set_dirs = abspath.(run_set_paths(name))
    if all(d -> d in set_dirs, abs_dirs)
      return name
    end
  end
  return nothing
end

function load_simulation_registry(path::AbstractString=SIMULATION_REGISTRY_PATH)
  isfile(path) || error("Simulation registry not found: $path")
  return TOML.parsefile(path)
end

function simulation_sets(registry=load_simulation_registry())
  return get(registry, "sets", Dict{String,Any}())
end

function list_run_sets(; registry=load_simulation_registry())
  return sort(collect(keys(simulation_sets(registry))))
end

function has_run_set(name::AbstractString; registry=load_simulation_registry())
  return haskey(simulation_sets(registry), String(name))
end

function get_run_set(name::AbstractString; registry=load_simulation_registry())
  sets = simulation_sets(registry)
  key = String(name)
  haskey(sets, key) || error("Unknown run set '$key'. Available: $(join(sort(collect(keys(sets))), ", "))")
  return sets[key]
end

function run_set_feedback(name::AbstractString; registry=load_simulation_registry())
  run_set = get_run_set(name; registry=registry)
  haskey(run_set, "thermodynamic_feedback") ||
    error("Run set '$name' does not declare thermodynamic_feedback")
  return Bool(run_set["thermodynamic_feedback"])
end

function project_path(path::AbstractString)
  isabspath(path) && return normpath(path)

  p = normpath(String(path))
  project_name = basename(SIMULATION_PROJECT_ROOT)
  p_unix = replace(p, '\\' => '/')

  # "Analysis/..." is already a workspace-root-relative path (the project name
  # is the first component); resolve it from the workspace root as-is.
  if p_unix == project_name || startswith(p_unix, project_name * "/")
    return normpath(joinpath(dirname(SIMULATION_PROJECT_ROOT), p))
  end

  # Otherwise the path is relative to the project root itself.
  return normpath(joinpath(SIMULATION_PROJECT_ROOT, p))
end

function run_set_entries(name::AbstractString; registry=load_simulation_registry())
  rs = get_run_set(name; registry=registry)
  return get(rs, "dirs", Any[])
end

function run_set_paths(name::AbstractString; registry=load_simulation_registry())
  return [project_path(String(entry["path"])) for entry in run_set_entries(name; registry=registry)]
end

function run_set_labels(name::AbstractString; registry=load_simulation_registry())
  return [String(get(entry, "label", basename(String(entry["path"])))) for entry in run_set_entries(name; registry=registry)]
end

function run_set_output_path(
  name::AbstractString,
  key::AbstractString;
  default::Union{Nothing,AbstractString}=nothing,
  registry=load_simulation_registry(),
)
  rs = get_run_set(name; registry=registry)
  outputs = get(rs, "outputs", Dict{String,Any}())
  if haskey(outputs, String(key))
    return project_path(String(outputs[String(key)]))
  end
  default === nothing && error("Run set '$name' has no output path '$key'")
  return project_path(default)
end

function read_run_config(run_dir::AbstractString)
  path = joinpath(run_dir, "config.toml")
  isfile(path) || error("Run config not found: $path")
  return TOML.parsefile(path)
end

function run_directory_feedback(run_dir::AbstractString)
  config = read_run_config(run_dir)
  haskey(config, "thermodynamic_feedback") || error(
    "Run config does not record thermodynamic_feedback: $run_dir. " *
    "Re-run with the current runner (run_anmbr_parallel.jl / run_ablation_theta.jl) so the mode is serialized."
  )
  return Bool(config["thermodynamic_feedback"])
end

function run_file_number(path::AbstractString)
  match_result = match(r"run(\d+)\.jld2$", basename(path))
  match_result === nothing && error("Cannot parse run number from: $path")
  return parse(Int, only(match_result.captures))
end

function _run_file_map(run_dir::AbstractString, subdir::AbstractString)
  root = joinpath(run_dir, subdir)
  isdir(root) || return Dict{Int,String}()
  files = filter(path -> endswith(lowercase(path), ".jld2"), readdir(root; join=true))
  mapping = Dict{Int,String}()
  for path in files
    run = run_file_number(path)
    haskey(mapping, run) && error("Duplicate run $run in $root")
    mapping[run] = path
  end
  return mapping
end

"""Return parameter/trajectory paths matched by the run number in each filename."""
function matched_run_files(run_dir::AbstractString; require_trajectories::Bool=true)
  run_dir = normpath(run_dir)
  param_map = _run_file_map(run_dir, "paras")
  isempty(param_map) && error("No parameter JLD2 files found in $(joinpath(run_dir, "paras"))")

  trajectory_map = _run_file_map(run_dir, "trajectories")
  param_runs = Set(keys(param_map))
  trajectory_runs = Set(keys(trajectory_map))
  if require_trajectories && param_runs != trajectory_runs
    missing = sort(collect(setdiff(param_runs, trajectory_runs)))
    extra = sort(collect(setdiff(trajectory_runs, param_runs)))
    error("Parameter/trajectory run mismatch in $run_dir; missing trajectories=$missing, extra trajectories=$extra")
  end

  runs = require_trajectories ? sort(collect(param_runs)) : sort(collect(union(param_runs, trajectory_runs)))
  return [(
    run=run,
    parameter_path=get(param_map, run, nothing),
    trajectory_path=get(trajectory_map, run, nothing),
  ) for run in runs]
end

function read_run_summary(run_dir::AbstractString)
  path = joinpath(run_dir, "summary.csv")
  isfile(path) || error("Run summary not found: $path")
  summary = CSV.read(path, DataFrame)
  required = Set([:run, :seed])
  missing = setdiff(required, Set(propertynames(summary)))
  isempty(missing) || error("Summary $path is missing columns: $(sort(collect(missing)))")
  return summary
end

function seed_run_map(run_dir::AbstractString)
  summary = read_run_summary(run_dir)
  mapping = Dict{Int,Int}()
  for row in eachrow(summary)
    seed = Int(row.seed)
    run = Int(row.run)
    haskey(mapping, seed) && error("Duplicate seed $seed in $(joinpath(run_dir, "summary.csv"))")
    mapping[seed] = run
  end
  return mapping
end

function run_files_for_seed(run_dir::AbstractString, seed::Integer)
  mapping = seed_run_map(run_dir)
  haskey(mapping, Int(seed)) || error("Seed $seed is not present in $run_dir")
  run = mapping[Int(seed)]
  files = Dict(item.run => item for item in matched_run_files(run_dir))
  haskey(files, run) || error("Summary maps seed $seed to run $run, but its JLD2 files are missing in $run_dir")
  return files[run]
end

function validate_seed_condition_grid(
  table::DataFrame,
  run_dirs::Vector{String};
  seed_column::Symbol=:seed,
  condition_column::Symbol=:SO4_mgL,
  config_key::String="sulfate_mg_L",
  label::String="Derived table",
)
  required = Set([seed_column, condition_column])
  missing = setdiff(required, Set(propertynames(table)))
  isempty(missing) || error("$label is missing columns: $(sort(collect(missing)))")

  expected = Set{Tuple{Int,Float64}}()
  for run_dir in run_dirs
    config = read_run_config(run_dir)
    haskey(config, config_key) || error("Run config has no '$config_key': $run_dir")
    condition = Float64(config[config_key])
    for seed in Int.(read_run_summary(run_dir).seed)
      push!(expected, (seed, condition))
    end
  end
  actual_rows = [
    (Int(row[seed_column]), Float64(row[condition_column]))
    for row in eachrow(table)
  ]
  actual = Set(actual_rows)
  length(actual_rows) == length(actual) || error("$label contains duplicate (seed, condition) rows")
  expected == actual || error(
    "$label does not match the selected raw run grid; " *
    "missing=$(length(setdiff(expected, actual))), extra=$(length(setdiff(actual, expected)))"
  )
  return true
end

function _same_struct_fields(left, right)
  typeof(left) == typeof(right) || return false
  return all(field -> getfield(left, field) == getfield(right, field), fieldnames(typeof(left)))
end

"""Validate that a seed represents the same reactions and strains in two conditions."""
function validate_paired_community(ps1, ps2; seed::Integer)
  (ps1.N, ps1.M, ps1.O) == (ps2.N, ps2.M, ps2.O) ||
    error("Seed $seed has incompatible model dimensions across conditions")
  for index in eachindex(ps1.reacs)
    _same_struct_fields(ps1.reacs[index], ps2.reacs[index]) ||
      error("Seed $seed reaction $index differs across conditions")
  end
  for index in eachindex(ps1.mics)
    _same_struct_fields(ps1.mics[index], ps2.mics[index]) ||
      error("Seed $seed strain $index differs across conditions")
  end
  return true
end

"""Validate that selected seeds retain identical communities across run directories."""
function validate_paired_run_directories(
  run_dirs::Vector{String};
  seeds::Union{Nothing,AbstractVector}=nothing,
)
  length(run_dirs) >= 2 || error("Paired-community validation requires at least two run directories")

  seed_maps = [seed_run_map(dir) for dir in run_dirs]
  seed_sets = [Set(keys(mapping)) for mapping in seed_maps]
  all(seed_set == first(seed_sets) for seed_set in seed_sets[2:end]) ||
    error("Paired run directories do not contain identical seed sets")

  selected_seeds = seeds === nothing ?
                   sort(collect(first(seed_sets))) :
                   sort(unique(Int.(seeds)))
  isempty(selected_seeds) && error("No seeds selected for paired-community validation")
  for seed in selected_seeds, (condition_index, seed_set) in enumerate(seed_sets)
    seed in seed_set || error("Seed $seed is absent from paired condition $condition_index")
  end

  parameter_paths = [
    Dict(
      item.run => item.parameter_path
      for item in matched_run_files(dir; require_trajectories=false)
      if item.parameter_path !== nothing
    )
    for dir in run_dirs
  ]
  for seed in selected_seeds
    reference_run = seed_maps[1][seed]
    reference_ps = load_parameters_silently(parameter_paths[1][reference_run])
    for condition_index in 2:length(run_dirs)
      run = seed_maps[condition_index][seed]
      ps = load_parameters_silently(parameter_paths[condition_index][run])
      validate_paired_community(reference_ps, ps; seed=seed)
    end
  end
  return true
end

function validate_run_directory(run_dir::AbstractString; require_trajectories::Bool=true)
  run_dir = normpath(run_dir)
  isdir(run_dir) || error("Simulation run directory does not exist: $run_dir")
  read_run_config(run_dir)
  summary = read_run_summary(run_dir)
  files = matched_run_files(run_dir; require_trajectories=require_trajectories)

  length(unique(Int.(summary.run))) == nrow(summary) ||
    error("summary.csv contains duplicate run IDs in $run_dir")
  length(unique(Int.(summary.seed))) == nrow(summary) ||
    error("summary.csv contains duplicate seeds in $run_dir")

  summary_runs = Set(Int.(summary.run))
  file_runs = Set(item.run for item in files if item.parameter_path !== nothing)
  summary_runs == file_runs || error(
    "summary.csv and parameter files disagree in $run_dir; " *
    "summary-only=$(sort(collect(setdiff(summary_runs, file_runs)))), " *
    "files-only=$(sort(collect(setdiff(file_runs, summary_runs))))"
  )

  if :solver_success in propertynames(summary) && any(.!Bool.(summary.solver_success))
    failed = Int.(summary.run[.!Bool.(summary.solver_success)])
    error("Run directory contains failed solver results: $run_dir, runs=$failed")
  end
  if :reached_tmax in propertynames(summary) && any(.!Bool.(summary.reached_tmax))
    incomplete = Int.(summary.run[.!Bool.(summary.reached_tmax)])
    error("Run directory contains trajectories that did not reach Tmax: $run_dir, runs=$incomplete")
  end

  config = read_run_config(run_dir)
  summary_by_run = Dict(Int(row.run) => Int(row.seed) for row in eachrow(summary))
  expected_feedback = haskey(config, "thermodynamic_feedback") ?
                      Bool(config["thermodynamic_feedback"]) :
                      nothing
  expected_schema = get(config, "schema_version", nothing)
  expected_revision = get(config, "model_revision", nothing)
  expected_hashes = get(config, "input_sha256", nothing)
  if haskey(config, "repeats")
    Int(config["repeats"]) == nrow(summary) ||
      error("config.toml repeats disagrees with summary.csv in $run_dir")
  end
  if haskey(config, "seed_sequence")
    order = sortperm(Int.(summary.run))
    actual_seed_sequence = Int.(summary.seed[order])
    expected_seed_sequence = Int.(config["seed_sequence"])
    actual_seed_sequence == expected_seed_sequence ||
      error("config.toml seed_sequence disagrees with summary.csv run order in $run_dir")
  end
  if expected_schema !== nothing && Int(expected_schema) >= 3
    expected_revision === nothing && error("Schema-$(expected_schema) config is missing model_revision: $run_dir")
    expected_hashes === nothing && error("Schema-$(expected_schema) config is missing input_sha256: $run_dir")
    expected_feedback === nothing && error("Schema-$(expected_schema) config is missing thermodynamic_feedback: $run_dir")
    haskey(config, "repeats") || error("Schema-$(expected_schema) config is missing repeats: $run_dir")
    haskey(config, "seed_sequence") || error("Schema-$(expected_schema) config is missing seed_sequence: $run_dir")
    haskey(config, "julia_version") || error("Schema-$(expected_schema) config is missing julia_version: $run_dir")
    haskey(config, "generated_at_utc") || error("Schema-$(expected_schema) config is missing generated_at_utc: $run_dir")
    sort(Int.(summary.run)) == collect(1:nrow(summary)) ||
      error("Schema-$(expected_schema) summary run IDs must be consecutive from 1 in $run_dir")

    required_summary = Set([:solver_retcode, :solver_success, :reached_tmax])
    missing_summary = setdiff(required_summary, Set(propertynames(summary)))
    isempty(missing_summary) || error(
      "Schema-$(expected_schema) summary is missing columns in $run_dir: $(sort(collect(missing_summary)))"
    )

    completion_path = joinpath(run_dir, "COMPLETED.toml")
    isfile(completion_path) || error("Schema-$(expected_schema) run is missing COMPLETED.toml: $run_dir")
    completion = TOML.parsefile(completion_path)
    required_completion = [
      "schema_version",
      "model_revision",
      "input_sha256",
      "completed_at_utc",
      "runs",
      "all_solver_success",
      "all_reached_tmax",
    ]
    missing_completion = filter(key -> !haskey(completion, key), required_completion)
    isempty(missing_completion) || error(
      "Schema-$(expected_schema) COMPLETED.toml is missing fields in $run_dir: $missing_completion"
    )
    Int(completion["schema_version"]) == Int(expected_schema) ||
      error("COMPLETED.toml schema_version disagrees with config.toml in $run_dir")
    String(completion["model_revision"]) == String(expected_revision) ||
      error("COMPLETED.toml model_revision disagrees with config.toml in $run_dir")
    completion["input_sha256"] == expected_hashes ||
      error("COMPLETED.toml input_sha256 disagrees with config.toml in $run_dir")
    Int(completion["runs"]) == nrow(summary) ||
      error("COMPLETED.toml run count disagrees with summary.csv in $run_dir")
    Bool(completion["all_solver_success"]) == all(Bool.(summary.solver_success)) ||
      error("COMPLETED.toml all_solver_success disagrees with summary.csv in $run_dir")
    Bool(completion["all_reached_tmax"]) == all(Bool.(summary.reached_tmax)) ||
      error("COMPLETED.toml all_reached_tmax disagrees with summary.csv in $run_dir")
    Bool(completion["all_solver_success"]) || error("COMPLETED.toml records solver failure in $run_dir")
    Bool(completion["all_reached_tmax"]) || error("COMPLETED.toml records an incomplete integration in $run_dir")
  end
  missing_metadata = Set{String}()

  function read_fields(path, names)
    return JLD2.jldopen(path, "r") do file
      Dict(name => (haskey(file, name) ? file[name] : nothing) for name in names)
    end
  end

  metadata_fields = ["seed", "schema_version", "model_revision", "input_sha256", "thermodynamic_feedback"]
  for item in files
    item.parameter_path === nothing && continue
    parameter = read_fields(item.parameter_path, metadata_fields)
    expected_seed = summary_by_run[item.run]
    if parameter["seed"] === nothing
      push!(missing_metadata, "parameter.seed")
    else
      Int(parameter["seed"]) == expected_seed ||
        error("Parameter seed mismatch in $(item.parameter_path): expected $expected_seed")
    end

    for field in ("schema_version", "model_revision", "input_sha256", "thermodynamic_feedback")
      parameter[field] === nothing && push!(missing_metadata, "parameter.$field")
    end
    parameter_schema = parameter["schema_version"] === nothing ? nothing : Int(parameter["schema_version"])
    if parameter_schema !== nothing && parameter_schema >= 3
      for field in metadata_fields
        parameter[field] === nothing &&
          error("Schema-$parameter_schema parameter file is missing '$field': $(item.parameter_path)")
      end
    end
    parameter["schema_version"] !== nothing && expected_schema !== nothing &&
      Int(parameter["schema_version"]) != Int(expected_schema) &&
      error("Parameter schema_version disagrees with config.toml in $(item.parameter_path)")
    parameter["model_revision"] !== nothing && expected_revision !== nothing &&
      String(parameter["model_revision"]) != String(expected_revision) &&
      error("Parameter model_revision disagrees with config.toml in $(item.parameter_path)")
    parameter["input_sha256"] !== nothing && expected_hashes !== nothing &&
      parameter["input_sha256"] != expected_hashes &&
      error("Parameter input_sha256 disagrees with config.toml in $(item.parameter_path)")
    parameter["thermodynamic_feedback"] !== nothing && expected_feedback !== nothing &&
      Bool(parameter["thermodynamic_feedback"]) != expected_feedback &&
      error("Parameter thermodynamic_feedback disagrees with config.toml in $(item.parameter_path)")

    item.trajectory_path === nothing && continue
    trajectory = read_fields(item.trajectory_path, metadata_fields)
    if trajectory["seed"] === nothing
      push!(missing_metadata, "trajectory.seed")
    else
      Int(trajectory["seed"]) == expected_seed ||
        error("Trajectory seed mismatch in $(item.trajectory_path): expected $expected_seed")
    end
    for field in ("schema_version", "model_revision", "input_sha256", "thermodynamic_feedback")
      trajectory[field] === nothing && push!(missing_metadata, "trajectory.$field")
      if parameter[field] !== nothing && trajectory[field] !== nothing
        parameter[field] == trajectory[field] ||
          error("Parameter/trajectory $field mismatch for run $(item.run) in $run_dir")
      end
    end
    if parameter_schema !== nothing && parameter_schema >= 3
      for field in metadata_fields
        trajectory[field] === nothing &&
          error("Schema-$parameter_schema trajectory file is missing '$field': $(item.trajectory_path)")
      end
    end

    if parameter_schema !== nothing
      ps = load_parameters_silently(item.parameter_path)
      terminal = read_fields(
        item.trajectory_path,
        ["final_pop", "final_conc", "final_energy", "final_phi"],
      )
      if parameter_schema >= 3
        for field in keys(terminal)
          terminal[field] === nothing &&
            error("Schema-$parameter_schema trajectory file is missing '$field': $(item.trajectory_path)")
        end
      end
      terminal["final_pop"] !== nothing && length(terminal["final_pop"]) != ps.N &&
        error("final_pop length disagrees with ps.N for run $(item.run) in $run_dir")
      terminal["final_conc"] !== nothing && length(terminal["final_conc"]) != ps.M &&
        error("final_conc length disagrees with ps.M for run $(item.run) in $run_dir")
      terminal["final_energy"] !== nothing && length(terminal["final_energy"]) != ps.N &&
        error("final_energy length disagrees with ps.N for run $(item.run) in $run_dir")
      terminal["final_phi"] !== nothing && length(terminal["final_phi"]) != ps.N &&
        error("final_phi length disagrees with ps.N for run $(item.run) in $run_dir")
    end
  end

  isempty(missing_metadata) || @warn(
    "Run directory has JLD2 files with missing semantic provenance",
    run_dir=run_dir,
    missing=sort(collect(missing_metadata)),
  )
  return true
end

function validate_run_set(
  name::AbstractString;
  require_trajectories::Bool=true,
  registry=load_simulation_registry(),
)
  run_set = get_run_set(name; registry=registry)
  status = String(get(run_set, "status", "unknown"))
  status == "broken" && error("Run set '$name' is marked broken: $(get(run_set, "notes", "no details"))")
  haskey(run_set, "thermodynamic_feedback") ||
    error("Run set '$name' must declare thermodynamic_feedback")

  entries = run_set_entries(name; registry=registry)
  isempty(entries) && error("Run set '$name' has no directories")
  reference_hashes = nothing
  reference_seeds = nothing
  reference_config = nothing
  paired_seeds = Bool(get(run_set, "paired_seeds", false))
  raw_varying_config_keys = get(run_set, "varying_config_keys", ["sulfate_mg_L"])
  raw_varying_config_keys isa AbstractVector ||
    error("Run set '$name' varying_config_keys must be an array")
  varying_config_keys = Set(String.(collect(raw_varying_config_keys)))
  union!(varying_config_keys, ["output_root", "generated_at_utc"])
  for (entry, path) in zip(entries, run_set_paths(name; registry=registry))
    validate_run_directory(path; require_trajectories=require_trajectories)
    config = read_run_config(path)
    if haskey(entry, "so4_mg_L") && haskey(config, "sulfate_mg_L")
      expected = Float64(entry["so4_mg_L"])
      actual = Float64(config["sulfate_mg_L"])
      isapprox(expected, actual; rtol=0.0, atol=1e-9) || error(
        "Run set '$name' labels $(basename(path)) as SO₄²⁻=$expected mg/L, but config.toml records $actual mg/L"
      )
    end
    if haskey(run_set, "thermodynamic_feedback")
      expected_feedback = Bool(run_set["thermodynamic_feedback"])
      if haskey(config, "thermodynamic_feedback")
        actual_feedback = Bool(config["thermodynamic_feedback"])
        expected_feedback == actual_feedback || error(
          "Run set '$name' feedback=$expected_feedback, but $(basename(path))/config.toml records $actual_feedback"
        )
      else
        error("Run set '$name' declares thermodynamic_feedback, but $(basename(path))/config.toml does not")
      end
    end
    if haskey(run_set, "schema_version")
      haskey(config, "schema_version") ||
        error("Run set '$name' declares a schema, but $(basename(path))/config.toml does not")
      Int(run_set["schema_version"]) == Int(config["schema_version"]) ||
        error("Run set '$name' schema_version disagrees with $(basename(path))/config.toml")
    end
    if haskey(run_set, "model_revision")
      haskey(config, "model_revision") ||
        error("Run set '$name' declares a model revision, but $(basename(path))/config.toml does not")
      String(run_set["model_revision"]) == String(config["model_revision"]) ||
        error("Run set '$name' model_revision disagrees with $(basename(path))/config.toml")
    end
    if haskey(config, "input_sha256")
      if reference_hashes === nothing
        reference_hashes = config["input_sha256"]
      else
        config["input_sha256"] == reference_hashes ||
          error("Run set '$name' mixes model/input SHA-256 maps")
      end
    end
    if paired_seeds
      seeds = Set(Int.(read_run_summary(path).seed))
      if reference_seeds === nothing
        reference_seeds = seeds
      else
        seeds == reference_seeds || error("Run set '$name' does not contain identical seed sets")
      end

      comparable_config = Dict(
        String(key) => value
        for (key, value) in config
        if !(String(key) in varying_config_keys)
      )
      if reference_config === nothing
        reference_config = comparable_config
      elseif comparable_config != reference_config
        all_keys = union(Set(keys(reference_config)), Set(keys(comparable_config)))
        differing_keys = sort([
          key for key in all_keys
                  if !haskey(reference_config, key) ||
                  !haskey(comparable_config, key) ||
                  reference_config[key] != comparable_config[key]
        ])
        error(
          "Run set '$name' changes undeclared config fields across paired conditions: " *
          "$(join(differing_keys, ", ")). Declare intentional variables in varying_config_keys."
        )
      end
    end
  end
  return true
end

function figure_output_dir(
  figure_name::AbstractString;
  run_set_name::Union{Nothing,AbstractString}=nothing,
  run_dirs::Vector{String}=String[],
  scope::Vector{String}=String[],
)
  source_id = if run_set_name !== nothing
    String(run_set_name)
  elseif isempty(run_dirs)
    "manual"
  else
    manual_source_id(run_dirs)
  end
  safe_scope = [replace(strip(part), r"[^A-Za-z0-9._=-]+" => "_") for part in scope]
  any(part -> isempty(part) || part in (".", ".."), safe_scope) &&
    error("Figure output scope contains an empty or unsafe path component: $scope")
  path = project_path(joinpath("results", source_id, String(figure_name), safe_scope...))
  mkpath(path)
  return path
end

function write_source_manifest(
  output_dir::AbstractString;
  script::AbstractString,
  run_set_name::Union{Nothing,AbstractString}=nothing,
  run_dirs::Vector{String}=String[],
  extra=Dict{String,Any}(),
  manifest_name::String="source.toml",
)
  all(key -> key isa AbstractString, keys(extra)) ||
    error("Source-manifest extra keys must be strings")
  basename(manifest_name) == manifest_name && endswith(manifest_name, ".toml") ||
    error("Source-manifest name must be a TOML filename without directories")
  reserved_keys = Set(["script", "script_path", "script_sha256", "generated_at_utc", "run_set", "run_dirs"])
  conflicting_keys = intersect(reserved_keys, Set(String.(collect(keys(extra)))))
  isempty(conflicting_keys) ||
    error("Source-manifest extra fields cannot replace reserved keys: $(sort(collect(conflicting_keys)))")
  manifest = Dict{String,Any}(
    "script" => String(script),
    "generated_at_utc" => string(now(UTC)),
    "run_set" => run_set_name === nothing ? "manual" : String(run_set_name),
    "run_dirs" => abspath.(run_dirs),
  )
  script_value = String(script)
  script_candidates = isabspath(script_value) ?
                      [normpath(script_value)] :
                      [
    project_path(script_value),
    joinpath(SIMULATION_PROJECT_ROOT, "src", "figures", basename(script_value)),
    joinpath(SIMULATION_PROJECT_ROOT, "src", "scripts", basename(script_value)),
  ]
  existing_scripts = unique(filter(isfile, script_candidates))
  if length(existing_scripts) == 1
    script_path = only(existing_scripts)
    manifest["script_path"] = replace(relpath(script_path, SIMULATION_PROJECT_ROOT), '\\' => '/')
    manifest["script_sha256"] = file_sha256(script_path)
  elseif length(existing_scripts) > 1
    error("Source script '$script' is ambiguous: $(join(existing_scripts, ", "))")
  end
  merge!(manifest, extra)
  mkpath(output_dir)
  open(joinpath(output_dir, manifest_name), "w") do io
    TOML.print(io, manifest)
  end
  return manifest
end

function available_run_set_names()
  return sort(collect(keys(simulation_sets())))
end

function resolve_run_set_args(
  args::Vector{String};
  default_set::Union{Nothing,String}=nothing,
  min_dirs::Int=1,
  require_trajectories::Bool=true,
)
  if isempty(args)
    default_set === nothing && error(
      "No run directories or --set <name> supplied. " *
      "Pick one of the registered run sets: $(join(available_run_set_names(), ", "))"
    )
    validate_run_set(default_set; require_trajectories=require_trajectories)
    return run_set_paths(default_set), default_set
  end

  if length(args) == 1 && has_run_set(args[1])
    validate_run_set(args[1]; require_trajectories=require_trajectories)
    return run_set_paths(args[1]), args[1]
  end

  if length(args) == 2 && args[1] in ("--set", "--run-set")
    validate_run_set(args[2]; require_trajectories=require_trajectories)
    return run_set_paths(args[2]), args[2]
  end

  if args[1] in ("--set", "--run-set")
    error("Expected exactly '--set <name>'; unexpected trailing arguments: $(join(args[3:end], " "))")
  end

  any(startswith(arg, "--") for arg in args) &&
    error("Unknown option in run-set arguments: $(join(filter(arg -> startswith(arg, "--"), args), ", "))")

  dirs = [project_path(arg) for arg in args]
  length(dirs) >= min_dirs || error("Expected at least $min_dirs run directories.")
  foreach(dir -> validate_run_directory(dir; require_trajectories=require_trajectories), dirs)
  return dirs, nothing
end

function resolve_single_run_args(
  args::Vector{String};
  default_set::Union{Nothing,String}=nothing,
  default_seed::Int=1,
)
  run_set_name = default_set
  condition = "1"
  seed = default_seed
  manual_dir = nothing
  explicit_set = false

  i = 1
  while i <= length(args)
    arg = args[i]
    if arg in ("--set", "--run-set")
      i += 1
      i <= length(args) || error("Expected a run-set name after $arg")
      run_set_name = args[i]
      explicit_set = true
    elseif arg == "--condition"
      i += 1
      i <= length(args) || error("Expected a label or index after --condition")
      condition = args[i]
    elseif arg == "--seed"
      i += 1
      i <= length(args) || error("Expected an integer after --seed")
      seed = parse(Int, args[i])
    elseif startswith(arg, "--")
      error("Unknown option: $arg")
    else
      manual_dir === nothing || error("Expected only one manual run directory")
      manual_dir = project_path(arg)
    end
    i += 1
  end

  if manual_dir !== nothing
    explicit_set && error("Do not combine a manual run directory with --set")
    validate_run_directory(manual_dir)
    return (run_dir=manual_dir, run_set_name=nothing, seed=seed, condition=basename(manual_dir))
  end

  run_set_name === nothing && error(
    "No run set supplied. Use --set <name> or pass a run directory. " *
    "Registered run sets: $(join(available_run_set_names(), ", "))"
  )
  validate_run_set(run_set_name)
  entries = run_set_entries(run_set_name)
  paths = run_set_paths(run_set_name)
  condition_index = tryparse(Int, condition)
  index = if condition_index !== nothing
    1 <= condition_index <= length(paths) || error("Condition index $condition_index is outside 1:$(length(paths))")
    condition_index
  else
    labels = run_set_labels(run_set_name)
    found = findfirst(==(condition), labels)
    found === nothing && error("Condition '$condition' not found. Available labels: $(join(labels, ", "))")
    found
  end
  label = String(get(entries[index], "label", basename(paths[index])))
  return (run_dir=paths[index], run_set_name=run_set_name, seed=seed, condition=label)
end

function print_run_set_summary(name::AbstractString)
  rs = get_run_set(name)
  println("Run set: ", name)
  println("  status: ", get(rs, "status", "unknown"))
  println("  ", get(rs, "description", ""))
  for (entry, path) in zip(run_set_entries(name), run_set_paths(name))
    label = get(entry, "label", basename(path))
    so4 = get(entry, "so4_mg_L", missing)
    so4_txt = ismissing(so4) ? "" : " SO₄²⁻=$(so4) mg/L"
    println("  - ", label, so4_txt, " -> ", path)
  end
end
