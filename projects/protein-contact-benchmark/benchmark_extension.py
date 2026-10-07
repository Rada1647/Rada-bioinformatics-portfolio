"""Execute a frozen multi-protein #332 benchmark without test-set tuning.

Inputs are prepared by prepare_extension.py. A validation selection is committed
to disk before any test arrays are loaded by the benchmark. Each trajectory is
processed separately. All inference intervals are conditional on a fixed clean
reference trajectory and artificial noise; they are not biological intervals.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import tempfile
import time

import numpy as np
import scipy
from scipy.stats import t as student_t

from methods import candidates, predict, predict_all_for_validation

ROOT = Path(__file__).resolve().parent
BASELINE_ORDER = ("raw", "linear", "hard_linear", "moving_average", "gaussian",
                  "median", "power_full", "power_matched")
DATASET_IDS = {"t4": 1, "ubiquitin": 2, "villin": 3}


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_new(path, value):
    """Never replace an existing output, including a completed scenario."""
    path = Path(path)
    encoded = json.dumps(value, indent=2, allow_nan=False) + "\n"
    with path.open("x", encoding="utf8") as stream:
        stream.write(encoded)


def key(params):
    # Protocol JSON may spell 1 while the implementation spells 1.0.
    # Numeric equality defines parameter identity, not its JSON spelling.
    return tuple(sorted(params.items()))


def seed_words(phase, dataset_id, run_index, scenario_index, noise_index):
    return [20261007, int(phase), int(dataset_id), int(run_index),
            int(scenario_index), int(noise_index)]


def contact_map(xyz, pairs, cutoff):
    """Float64 distances and strict threshold, with bounded pair temporaries."""
    answer = np.empty((len(xyz), len(pairs)), dtype=bool)
    for start in range(0, len(pairs), 256):
        block = pairs[start:start + 256]
        delta = (xyz[:, block[:, 0]].astype(np.float64)
                 - xyz[:, block[:, 1]].astype(np.float64))
        answer[:, start:start + 256] = np.linalg.norm(delta, axis=2) < cutoff
    return answer


def corrupt(X, xyz, pairs, scenario, words):
    rng = np.random.default_rng(np.random.SeedSequence(words))
    kind, epsilon = scenario["kind"], float(scenario["level"])
    if kind == "bitflip":
        return np.logical_xor(X, rng.random(X.shape) < epsilon)
    if kind == "correlated":
        rho = float(scenario["rho"])
        errors = np.empty(X.shape, dtype=bool)
        errors[0] = rng.random(X.shape[1]) < epsilon
        p01, p11 = epsilon * (1 - rho), 1 - (1 - epsilon) * (1 - rho)
        for i in range(1, len(X)):
            errors[i] = rng.random(X.shape[1]) < np.where(errors[i - 1], p11, p01)
        return np.logical_xor(X, errors)
    if kind == "jitter":
        if xyz is None:
            raise ValueError("coordinate jitter requires CA coordinates")
        jittered = xyz.astype(np.float64) + rng.normal(0, epsilon, size=xyz.shape)
        return contact_map(jittered, pairs, float(scenario["cutoff_A"]))
    raise ValueError(f"unknown corruption kind: {kind}")


def strata(X, frequency):
    truth = np.asarray(X, dtype=bool)
    rare = frequency <= 0.10
    transition = np.zeros_like(truth)
    changes = truth[1:] != truth[:-1]
    transition[1:] |= changes
    transition[:-1] |= changes
    brief = np.zeros_like(truth)
    brief_event_count = 0
    for j in range(X.shape[1]):
        edges = np.diff(np.r_[False, truth[:, j], False].astype(np.int8))
        starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
        for first, stop in zip(starts, ends):
            if stop - first <= 3 and first > 0 and stop < len(X):
                brief[first:stop, j] = True
                brief_event_count += 1
    return {"rare_features": rare, "rare_positive": truth & rare[None, :],
            "transition": transition, "short_positive": brief,
            "brief_event_count": brief_event_count}


def stratum_counts(masks):
    return {name: int(mask) if name == "brief_event_count" else int(mask.sum())
            for name, mask in masks.items()}


def selected_mean(array, mask):
    return float(array[mask].mean()) if np.any(mask) else None


def metrics(P, X, masks):
    err = P - X
    square = err * err
    binary, truth = P >= 0.5, X.astype(bool)
    tp = int(np.count_nonzero(binary & truth))
    fp = int(np.count_nonzero(binary & ~truth))
    fn = int(np.count_nonzero(~binary & truth))
    return {
        "brier": float(square.mean()), "mae": float(np.abs(err).mean()),
        "occupancy_mae": float(np.abs(P.mean(axis=0) - X.mean(axis=0)).mean()),
        "occupancy_signed_bias": float(err.mean()),
        "precision": tp / (tp + fp) if tp + fp else None,
        "recall": tp / (tp + fn) if tp + fn else None,
        "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None,
        "rare_brier": (float(square[:, masks["rare_features"]].mean())
                       if masks["rare_features"].any() else None),
        "rare_positive_recall": selected_mean(binary, masks["rare_positive"]),
        "short_positive_recall": selected_mean(binary, masks["short_positive"]),
        "transition_brier": selected_mean(square, masks["transition"]),
    }


def available_mean(values):
    """Protocol v1.1: describe available strata and separately report counts."""
    available = [value for value in values if value is not None]
    return float(np.mean(available)) if available else None


def mean_metrics(metric_rows):
    return {name: available_mean([row[name] for row in metric_rows])
            for name in metric_rows[0]}


def metric_counts(metric_rows):
    return {name: sum(row[name] is not None for row in metric_rows)
            for name in metric_rows[0]}


def paired_interval(rows, baseline):
    a = np.asarray([row["methods"][baseline]["brier"] for row in rows])
    b = np.asarray([row["methods"]["paper332"]["brier"] for row in rows])
    difference = a - b
    half = (float(student_t.ppf(0.975, len(difference) - 1)
                  * difference.std(ddof=1) / np.sqrt(len(difference)))
            if len(difference) > 1 else None)
    mean = float(difference.mean())
    return {
        "baseline": baseline, "n_noise_realizations": len(difference),
        "baseline_mean_brier": float(a.mean()), "paper332_mean_brier": float(b.mean()),
        "mean_brier_gain": mean,
        "relative_mean_brier_gain": mean / float(a.mean()) if a.mean() > 0 else None,
        "noise_only_paired_95pct_ci": [mean - half, mean + half] if half is not None else None,
        "positive_noise_differences": int(np.count_nonzero(difference > 0)),
        "negative_noise_differences": int(np.count_nonzero(difference < 0)),
        "interpretation": "Pointwise descriptive interval conditional on this fixed trajectory; not biological uncertainty; no multiplicity-adjusted significance claim.",
    }


def composite_success(per_run, aggregate, baseline, rules):
    baseline_brier = aggregate[baseline]["brier"]
    gain = baseline_brier - aggregate["paper332"]["brier"]
    relative = gain / baseline_brier if baseline_brier > 0 else None
    run_conditions = []
    for run in per_run:
        b, p = run["aggregate"][baseline], run["aggregate"]["paper332"]
        rare = (b["rare_positive_recall"] - p["rare_positive_recall"]
                if b["rare_positive_recall"] is not None and p["rare_positive_recall"] is not None else None)
        brief = (b["short_positive_recall"] - p["short_positive_recall"]
                 if b["short_positive_recall"] is not None and p["short_positive_recall"] is not None else None)
        run_conditions.append({
            "run_id": run["run_id"], "brier_gain": b["brier"] - p["brier"],
            "positive_brier_gain": bool(b["brier"] > p["brier"]),
            "rare_recall_loss": rare, "brief_recall_loss": brief,
            "rare_guardrail_pass": None if rare is None else bool(rare <= rules["max_rare_recall_loss_each_run"]),
            "brief_guardrail_pass": None if brief is None else bool(brief <= rules["max_brief_recall_loss_each_run"]),
        })
    conditions = {
        "relative_gain_at_least_minimum": None if relative is None else bool(relative >= rules["relative_mean_brier_gain_min"]),
        "positive_gain_every_test_run": all(row["positive_brier_gain"] for row in run_conditions),
        "rare_guardrail_every_test_run": (None if any(row["rare_guardrail_pass"] is None for row in run_conditions)
                                           else all(row["rare_guardrail_pass"] for row in run_conditions)),
        "brief_guardrail_every_test_run": (None if any(row["brief_guardrail_pass"] is None for row in run_conditions)
                                            else all(row["brief_guardrail_pass"] for row in run_conditions)),
    }
    status = ("not_evaluable" if any(v is None for v in conditions.values())
              else "pass" if all(conditions.values()) else "fail")
    return {"baseline": baseline, "status": status, "relative_mean_brier_gain": relative,
            "conditions": conditions, "per_run": run_conditions,
            "scope": "Prespecified practical composite for the primary scenario only; other scenarios descriptive."}


def resolve_input(value, manifest_path):
    path = Path(value)
    if path.is_absolute():
        return path
    for candidate in (manifest_path.parent / path, ROOT / path, Path.cwd() / path):
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(f"cannot resolve input {value!r} from {manifest_path}")


def load_run(run, cutoff, manifest_path, need_xyz):
    path = resolve_input(run["contacts"][str(int(cutoff))], manifest_path)
    with np.load(path, allow_pickle=False) as data:
        X, pairs = data["X"], data["pairs"]
        frequency = data["frequency"] if "frequency" in data else data["discovery_frequency"]
        times = data["time_ps"] if "time_ps" in data else None
    if X.ndim != 2 or len(X) < 1 or X.shape[1] < 1:
        raise ValueError(f"non-evaluable empty contact array: {path}")
    if not np.all((X == 0) | (X == 1)):
        raise ValueError(f"nonbinary contact array: {path}")
    if pairs.shape != (X.shape[1], 2) or frequency.shape != (X.shape[1],):
        raise ValueError(f"contact metadata mismatch: {path}")
    if not np.issubdtype(pairs.dtype, np.integer) or not np.all(pairs[:, 1] - pairs[:, 0] >= 4):
        raise ValueError(f"invalid nonlocal CA pairs: {path}")
    if not np.all(np.isfinite(frequency)) or not np.all((frequency > 0) & (frequency < 1)):
        raise ValueError(f"discovery-constant feature was selected: {path}")
    if times is not None and (times.shape != (len(X),) or not np.all(np.isfinite(times))
                              or not np.all(np.diff(times) > 0)):
        raise ValueError(f"invalid chronological times: {path}")
    xyz = None
    if need_xyz:
        with np.load(resolve_input(run["coordinates_path"], manifest_path), allow_pickle=False) as data:
            xyz = data["xyz"]
        if xyz.ndim != 3 or xyz.shape[0] != len(X) or xyz.shape[2] != 3 or not np.all(np.isfinite(xyz)):
            raise ValueError(f"invalid coordinate array for {run['id']}")
        if pairs.min() < 0 or pairs.max() >= xyz.shape[1]:
            raise ValueError(f"pair atom index is outside coordinate array for {run['id']}")
    return X.astype(bool), pairs, frequency.astype(float), xyz, times


def validate_compatible(reference, pairs, frequency):
    if reference is None:
        return pairs.copy(), frequency.copy()
    if not np.array_equal(reference[0], pairs) or not np.array_equal(reference[1], frequency):
        raise ValueError("selected features or discovery frequencies differ between runs")
    return reference


def feature_identity(pairs, frequency):
    return {"pairs_sha256": hashlib.sha256(np.asarray(pairs, dtype="<i8").tobytes()).hexdigest(),
            "discovery_frequency_sha256": hashlib.sha256(np.asarray(frequency, dtype="<f8").tobytes()).hexdigest()}


def eligible_grids(protocol, scenario):
    grids = protocol["methods"]
    available = candidates()
    if grids != available:
        raise ValueError("methods.py candidate grids do not exactly match the supplied protocol")
    return {name: values for name, values in grids.items()
            if name != "oracle_debiased_linear" or scenario["kind"] in ("bitflip", "correlated")}


def oracle_epsilon(scenario):
    return scenario["level"] if scenario["kind"] in ("bitflip", "correlated") else None


def common_identity(protocol_path, manifest_path, dataset, scenario, scenario_index):
    return {"dataset": dataset, "scenario": scenario, "scenario_index": scenario_index,
            "protocol_sha256": digest(protocol_path), "prepared_manifest_sha256": digest(manifest_path),
            "code_sha256": {"benchmark_extension.py": digest(__file__), "methods.py": digest(ROOT / "methods.py")}}


def validate_selection_identity(selection, identity):
    for name, value in identity.items():
        if selection.get(name) != value:
            raise ValueError(f"saved validation identity mismatch for {name}; use a new output directory")


def validation_phase(protocol, manifest, manifest_path, identity, out):
    start = time.perf_counter()
    scenario = identity["scenario"]
    path = out / f"{scenario['name']}_validation.json"
    if path.exists():
        raise FileExistsError(f"will not overwrite saved validation: {path}")
    grids = eligible_grids(protocol, scenario)
    lookup = {name: {key(params): i for i, params in enumerate(params_list)}
              for name, params_list in grids.items()}
    sums = {name: np.zeros(len(params_list)) for name, params_list in grids.items()}
    records, reference = [], None
    runs = [run for run in manifest["runs"] if run["split"] == "validation"]
    if not runs:
        raise ValueError("no validation trajectories")
    for run in runs:
        X, pairs, freq, xyz, _ = load_run(run, scenario["cutoff_A"], manifest_path, scenario["kind"] == "jitter")
        reference = validate_compatible(reference, pairs, freq)
        for index in range(protocol["noise"]["validation_count"]):
            words = seed_words(1, manifest["dataset_id"], run["run_index"], identity["scenario_index"], index)
            Z = corrupt(X, xyz, pairs, scenario, words)
            scores = {name: [None] * len(params_list) for name, params_list in grids.items()}
            for name, params, P in predict_all_for_validation(Z, oracle_epsilon(scenario)):
                slot = lookup[name][key(params)]
                score = float(np.square(P - X).mean())
                scores[name][slot] = score
                sums[name][slot] += score
            if any(any(value is None for value in values) for values in scores.values()):
                raise RuntimeError("validation generator omitted a declared candidate")
            records.append({"run_id": run["id"], "run_index": run["run_index"],
                            "noise_index": index, "seed_words": words, "candidate_brier": scores})
        print(f"validation {identity['dataset']} {scenario['name']} run={run['id']} complete", flush=True)
    means = {name: (values / len(records)).tolist() for name, values in sums.items()}
    settings = {name: params_list[int(np.argmin(means[name]))] for name, params_list in grids.items()}
    best_scores = {name: min(values) for name, values in means.items()}
    baseline = min(BASELINE_ORDER, key=lambda name: best_scores[name])
    matched = min((name for name in BASELINE_ORDER if name != "power_full"), key=lambda name: best_scores[name])
    result = {**identity, "status": "complete", "completed_utc": utcnow(),
              "feature_identity": feature_identity(*reference),
              "candidate_grids": grids, "candidate_counts": {name: len(values) for name, values in grids.items()},
              "records": records, "mean_candidate_brier": means, "selected_parameters": settings,
              "selected_validation_brier": best_scores, "selected_baseline": baseline,
              "selected_matched_baseline": matched,
              "aggregation": "Equal weight to each trajectory and each of its identically many noise realizations.",
              "seconds": time.perf_counter() - start}
    save_new(path, result)
    return result


def test_phase(protocol, manifest, manifest_path, identity, selection, out):
    start = time.perf_counter()
    scenario = identity["scenario"]
    path = out / f"{scenario['name']}_test.json"
    if path.exists():
        raise FileExistsError(f"will not overwrite completed tests: {path}")
    validate_selection_identity(selection, identity)
    selected = selection["selected_parameters"]
    methods = list(selected)
    per_run, reference = [], None
    method_seconds = {name: 0.0 for name in methods}
    runs = [run for run in manifest["runs"] if run["split"] == "test"]
    if not runs:
        raise ValueError("no test trajectories")
    for run in runs:
        X, pairs, freq, xyz, times = load_run(run, scenario["cutoff_A"], manifest_path, scenario["kind"] == "jitter")
        reference = validate_compatible(reference, pairs, freq)
        if feature_identity(pairs, freq) != selection["feature_identity"]:
            raise ValueError("test feature identity differs from saved validation")
        masks, rows = strata(X, freq), []
        for index in range(protocol["noise"]["test_count"]):
            words = seed_words(2, manifest["dataset_id"], run["run_index"], identity["scenario_index"], index)
            Z = corrupt(X, xyz, pairs, scenario, words)
            result = {"noise_index": index, "seed_words": words, "methods": {}}
            for name in methods:
                tic = time.perf_counter()
                P = predict(Z, name, selected[name], noise_level=oracle_epsilon(scenario))
                result["methods"][name] = metrics(P, X, masks)
                if name == "paper332":
                    cubic_binary = P >= 0.5
                del P
                method_seconds[name] += time.perf_counter() - tic
            H = predict(Z, "linear", {"t": selected["paper332"]["t"]})
            result["same_t_linear"] = metrics(H, X, masks)
            result["same_t_threshold_disagreements"] = int(np.count_nonzero((H >= 0.5) != cubic_binary))
            del H, cubic_binary, Z
            rows.append(result)
        run_result = {
            "run_id": run["id"], "run_index": run["run_index"], "shape": list(X.shape),
            "sampled_time_range_ps": None if times is None else [float(times[0]), float(times[-1])],
            "stratum_counts": stratum_counts(masks), "noise_rows": rows,
            "aggregate": {name: mean_metrics([row["methods"][name] for row in rows]) for name in methods},
            "aggregate_contributing_noise_count": {name: metric_counts([row["methods"][name] for row in rows]) for name in methods},
            "same_t_linear_aggregate": mean_metrics([row["same_t_linear"] for row in rows]),
            "same_t_linear_contributing_noise_count": metric_counts([row["same_t_linear"] for row in rows]),
            "same_t_threshold_disagreements": sum(row["same_t_threshold_disagreements"] for row in rows),
            "paired_comparisons": {name: paired_interval(rows, name) for name in methods if name != "paper332"},
        }
        per_run.append(run_result)
        print(f"test {identity['dataset']} {scenario['name']} run={run['id']} complete", flush=True)
    aggregate = {name: mean_metrics([run["aggregate"][name] for run in per_run]) for name in methods}
    result = {**identity, "status": "complete", "completed_utc": utcnow(),
              "validation_file_sha256": digest(out / f"{scenario['name']}_validation.json"),
              "dataset_role": protocol["datasets"][identity["dataset"]]["role"],
              "selected_parameters": selected, "selected_baseline": selection["selected_baseline"],
              "selected_matched_baseline": selection["selected_matched_baseline"],
              "per_run": per_run, "aggregate": aggregate,
              "aggregate_contributing_run_count": {name: metric_counts([run["aggregate"][name] for run in per_run]) for name in methods},
              "same_t_linear_aggregate": mean_metrics([run["same_t_linear_aggregate"] for run in per_run]),
              "same_t_linear_contributing_run_count": metric_counts([run["same_t_linear_aggregate"] for run in per_run]),
              "same_t_threshold_disagreements": sum(run["same_t_threshold_disagreements"] for run in per_run),
              "method_prediction_and_metrics_seconds": method_seconds,
              "seconds": time.perf_counter() - start,
              "aggregation": "Available noise means within each trajectory, followed by equally weighted available trajectory means. Null strata are retained per run; contributing counts are explicit. Missing guardrail strata make practical success not fully evaluable.",
              "interpretation": protocol["uncertainty"]}
    if scenario.get("primary"):
        criterion_key = "stress_success_descriptive" if identity["dataset"] == "villin" else "primary_success"
        result[criterion_key] = composite_success(per_run, aggregate, selection["selected_baseline"], protocol["primary_success"])
        result["matched_budget_success_secondary"] = composite_success(per_run, aggregate, selection["selected_matched_baseline"], protocol["primary_success"])
        if identity["dataset"] == "villin":
            for criterion in (result[criterion_key], result["matched_budget_success_secondary"]):
                criterion["scope"] = "Descriptive adaptive-sampling stress criterion only; not primary confirmation or biological replication."
    save_new(path, result)
    return result


def clean_phase(protocol, manifest, manifest_path, identity, selection, out):
    path = out / "primary_clean.json"
    if path.exists():
        raise FileExistsError(f"will not overwrite clean-control output: {path}")
    start = time.perf_counter()
    validate_selection_identity(selection, identity)
    rows, reference = [], None
    for run in manifest["runs"]:
        if run["split"] != "test":
            continue
        X, pairs, freq, _, _ = load_run(run, identity["scenario"]["cutoff_A"], manifest_path, False)
        reference = validate_compatible(reference, pairs, freq)
        if feature_identity(pairs, freq) != selection["feature_identity"]:
            raise ValueError("clean-control feature identity differs from saved validation")
        masks = strata(X, freq)
        row = {"run_id": run["id"], "run_index": run["run_index"], "shape": list(X.shape),
               "stratum_counts": stratum_counts(masks), "methods": {}}
        for name, params in selection["selected_parameters"].items():
            P = predict(X, name, params, noise_level=identity["scenario"]["level"])
            row["methods"][name] = metrics(P, X, masks)
            del P
        rows.append(row)
    result = {**identity, "status": "complete", "completed_utc": utcnow(),
              "selected_parameters": selection["selected_parameters"], "per_run": rows,
              "aggregate": {name: mean_metrics([row["methods"][name] for row in rows])
                            for name in selection["selected_parameters"]},
              "aggregate_contributing_run_count": {name: metric_counts([row["methods"][name] for row in rows])
                                                   for name in selection["selected_parameters"]},
              "oracle_epsilon": identity["scenario"]["level"],
              "oracle_caveat": "Uses epsilon from noisy primary validation on clean input, intentionally exposing misspecification; no clean-input retuning.",
              "seconds": time.perf_counter() - start}
    save_new(path, result)
    return result


def run_benchmark(args):
    protocol_path, manifest_path = args.protocol.resolve(), args.manifest.resolve()
    protocol = json.loads(protocol_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("protocol_sha256") != digest(protocol_path):
        raise ValueError("prepared manifest does not match the supplied frozen protocol")
    if manifest["dataset_id"] != DATASET_IDS[args.dataset]:
        raise ValueError("dataset name and dataset_id disagree")
    ids = [run["id"] for run in manifest["runs"]]
    indices = [run["run_index"] for run in manifest["runs"]]
    if len(ids) != len(set(ids)) or len(indices) != len(set(indices)):
        raise ValueError("run IDs and seed run indices must be unique")
    scenarios = [(index, scenario) for index, scenario in enumerate(protocol["scenarios"])
                 if args.scenario in ("all", scenario["name"])]
    if not scenarios:
        raise ValueError(f"unknown scenario {args.scenario}")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    # Reject completed outputs before doing any scenario work in this invocation.
    for _, scenario in scenarios:
        suffixes = ("validation", "test") if args.phase == "all" else (args.phase,)
        for suffix in suffixes:
            target = out / f"{scenario['name']}_{suffix}.json"
            if target.exists():
                raise FileExistsError(f"will not overwrite {target}")
    files = {manifest_path, protocol_path, Path(__file__).resolve(), ROOT / "methods.py"}
    for run in manifest["runs"]:
        for value in [run["coordinates_path"], *run["contacts"].values()]:
            files.add(resolve_input(value, manifest_path))
    invocation = {"started_utc": utcnow(), "dataset": args.dataset, "phase": args.phase,
                  "scenario": args.scenario, "command": sys.argv,
                  "inputs_and_code_sha256": {str(path): digest(path) for path in sorted(files)},
                  "environment": {"python": platform.python_version(), "numpy": np.__version__,
                                  "scipy": scipy.__version__, "platform": platform.platform()}}
    # A manifest hash alone does not detect a subsequently edited input file.
    # Check the prepared manifest's individual data hashes before execution.
    prepared_hashes = {}
    for run in manifest["runs"]:
        coords = resolve_input(run["coordinates_path"], manifest_path)
        actual_coords = invocation["inputs_and_code_sha256"][str(coords)]
        if run.get("coordinates_sha256", actual_coords) != actual_coords:
            raise ValueError(f"prepared coordinate checksum mismatch: {run['id']}")
        hashes = {"coordinates": actual_coords, "contacts": {}}
        for cutoff, value in run["contacts"].items():
            actual = invocation["inputs_and_code_sha256"][str(resolve_input(value, manifest_path))]
            if run.get("contacts_sha256", {}).get(cutoff, actual) != actual:
                raise ValueError(f"prepared contact checksum mismatch: {run['id']} cutoff {cutoff}")
            hashes["contacts"][cutoff] = actual
        prepared_hashes[run["id"]] = hashes
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    save_new(out / f"invocation_{timestamp}.json", invocation)
    for scenario_index, scenario in scenarios:
        identity = common_identity(protocol_path, manifest_path, args.dataset, scenario, scenario_index)
        identity["prepared_inputs_sha256"] = prepared_hashes
        if args.phase in ("validation", "all"):
            selection = validation_phase(protocol, manifest, manifest_path, identity, out)
        else:
            selection = json.loads((out / f"{scenario['name']}_validation.json").read_text())
            validate_selection_identity(selection, identity)
        if args.phase in ("test", "all"):
            result = test_phase(protocol, manifest, manifest_path, identity, selection, out)
            if scenario.get("primary"):
                clean_phase(protocol, manifest, manifest_path, identity, selection, out)
            print(json.dumps({"dataset": args.dataset, "scenario": scenario["name"],
                              "brier": {name: row["brier"] for name, row in result["aggregate"].items()},
                              "baseline": selection["selected_baseline"],
                              "composite": result.get("primary_success", result.get("stress_success_descriptive", {})).get("status"),
                              "dataset_role": result["dataset_role"]}), flush=True)


def self_test():
    """Synthetic checks only; does not inspect or evaluate protein data."""
    X = np.zeros((14, 3), dtype=bool)
    X[0:2, 0] = True                  # boundary-censored, must be excluded
    X[4:7, 0] = True                  # complete three-frame event
    X[8:12, 1] = True                 # complete four-frame event, excluded
    X[12:, 2] = True                  # other boundary-censored event
    masks = strata(X, np.asarray([0.08, 0.5, 0.01]))
    assert masks["short_positive"].sum() == 3 and masks["brief_event_count"] == 1
    assert metrics(X.astype(float), X, masks)["brier"] == 0
    assert metrics(X.astype(float), X, masks)["short_positive_recall"] == 1
    empty = np.zeros((5, 1), dtype=bool)
    em = metrics(empty.astype(float), empty, strata(empty, np.asarray([0.4])))
    assert em["rare_positive_recall"] is None and em["short_positive_recall"] is None
    assert available_mean([1.0, None]) == 1.0 and available_mean([None]) is None
    assert metric_counts([{"a": None}, {"a": 1.0}]) == {"a": 1}
    guardrail_rules = {"relative_mean_brier_gain_min": 0.02,
                       "max_rare_recall_loss_each_run": 0.02,
                       "max_brief_recall_loss_each_run": 0.02}
    controls = {"raw": {"brier": 0.1, "rare_positive_recall": 0.8, "short_positive_recall": 0.7},
                "paper332": {"brier": 0.09, "rare_positive_recall": 0.81, "short_positive_recall": 0.71}}
    fake_runs = [{"run_id": "synthetic", "aggregate": controls}]
    assert composite_success(fake_runs, controls, "raw", guardrail_rules)["status"] == "pass"
    controls["paper332"]["rare_positive_recall"] = 0.7
    assert composite_success(fake_runs, controls, "raw", guardrail_rules)["status"] == "fail"
    controls["paper332"]["short_positive_recall"] = None
    assert composite_success(fake_runs, controls, "raw", guardrail_rules)["status"] == "not_evaluable"
    scenario = {"kind": "correlated", "level": 0.1, "rho": 0.8, "cutoff_A": 8}
    errors = corrupt(np.zeros((50000, 12), bool), None, None, scenario, [1, 2, 3])
    prevalence = errors.mean()
    p01 = errors[1:][~errors[:-1]].mean()
    p10 = (~errors[1:][errors[:-1]]).mean()
    assert abs(prevalence - 0.1) < 0.005, prevalence
    assert abs(p01 - 0.02) < 0.002 and abs(p10 - 0.18) < 0.006, (p01, p10)
    assert np.array_equal(corrupt(X, None, None, scenario, [4]), corrupt(X, None, None, scenario, [4]))
    # A synthetic five-run end-to-end check uses the actual frozen grids/seeds.
    with tempfile.TemporaryDirectory(prefix="contact332_synthetic_") as tmp:
        tmp = Path(tmp)
        protocol_path = ROOT / "protocol" / "protocol_v1.1.json"
        protocol = json.loads(protocol_path.read_text())
        runs = []
        pairs = np.asarray([[0, 4], [0, 5], [1, 5]], dtype=int)
        freq = np.asarray([0.08, 0.5, 0.5])
        for i, split in enumerate(("discovery", "validation", "test", "test", "test")):
            rng = np.random.default_rng(400 + i)
            xyz = rng.normal(0, 3, size=(24, 6, 3)).astype(np.float32)
            coords = tmp / f"run{i}_coords.npz"
            np.savez(coords, xyz=xyz)
            paths = {}
            for cutoff in (7, 8, 9):
                contact = tmp / f"run{i}_contact{cutoff}.npz"
                np.savez(contact, X=contact_map(xyz, pairs, cutoff), pairs=pairs,
                         frequency=freq, time_ps=np.arange(24) * 200)
                paths[str(cutoff)] = str(contact)
            runs.append({"id": f"synthetic{i}", "run_index": i, "split": split,
                         "coordinates_path": str(coords), "contacts": paths})
        manifest_path = tmp / "manifest.json"
        save_new(manifest_path, {"dataset_id": 1, "protocol_sha256": digest(protocol_path), "runs": runs})
        args = argparse.Namespace(dataset="t4", protocol=protocol_path, manifest=manifest_path,
                                  scenario="all", phase="all", out=tmp / "results")
        run_benchmark(args)
        for scenario in protocol["scenarios"]:
            result = json.loads((args.out / f"{scenario['name']}_test.json").read_text())
            assert len(result["per_run"]) == 3
            assert all(len(run["noise_rows"]) == 16 for run in result["per_run"])
            assert result["same_t_threshold_disagreements"] >= 0
        try:
            run_benchmark(args)
        except FileExistsError:
            pass
        else:
            raise AssertionError("completed-output overwrite was not rejected")
    print(json.dumps({"synthetic_self_test": "pass", "correlated_prevalence": float(prevalence),
                      "correlated_p01": float(p01), "correlated_p10": float(p10)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=tuple(DATASET_IDS))
    parser.add_argument("--protocol", type=Path, default=ROOT / "protocol" / "protocol_v1.1.json")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--scenario", default="all")
    parser.add_argument("--phase", choices=("validation", "test", "all"), default="all")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.dataset is None:
        parser.error("--dataset is required unless --self-test is used")
    if args.manifest is None:
        args.manifest = ROOT / "data" / "prepared" / args.dataset / "manifest.json"
    if args.out is None:
        args.out = ROOT / "results" / args.dataset
    run_benchmark(args)


if __name__ == "__main__":
    main()
