# Project Instructions

For AI coding agents working in this repository:

- After making code changes, commit them to git before finishing the task.
- Push the commit to the configured remote branch when the push is expected to work.
- Do not commit unrelated user changes. Stage only files changed for the current task.
- If a push fails because credentials, network access, or permissions are unavailable, report the commit hash and the exact push failure.
- Use concise commit messages that describe the behavior changed, not just the files edited.

## Project Notes

- `predictionBreak.m` is the main analysis class for this repo. Prefer scoped edits there and in the specific figure script the user names.
- Spike-rate bins use `options.binw` in milliseconds. Convert counts/bin to Hz with `counts / (options.binw/1000)`, not `counts / (options.binw*1000)`.
- Spike analyses are cached as `spikes_*.mat` in `obj.cacheDir`. If rate units, binning, tuning filters, or event alignment change, tell the user to delete or regenerate those caches with `is_cache=false`.
- `res` from `analyze_spikes` is a 7-row cell array: cluster id, rate matrix, time vector, raw event-aligned spike times, p-value, cluster group, template.
- For `flip_time`, spike tuning p-values are intentionally replaced with `escape_time` p-values via `get_spikes_tuned_units_for_flip_trials`, because flip trials can be sparse.
- `plot_spikes_raster_plot` accepts a single recording name or a cell array of recordings. When `cluster_id=nan`, it selects tuned units (`p < alpha`) if `is_tuned=true`; otherwise it plots all clusters.
- `plot_avg_spikes_rate` supports `smooth_kernel`. The smoothing is causal half-Gaussian, not symmetric Savitzky-Golay, to avoid pre-event leakage.
- Be careful with file handles. Prefer `currentDataObj.closeOpenFiles()` through `close_current_data_files`; avoid broad `fclose('all')` inside normal per-record loops unless there is a clear reason.
- `figure6.m` writes per-observation statistics to `figure6_stats.csv` in `pb.outputDir`, with columns `animal_id`, `metric`, `label`, and `power_diff`.
