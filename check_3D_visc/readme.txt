for run_process_visc_all_layers.sh, the output files' layer id is counted from bottom (1) to top (global layers number), i.e., same order as that used in SVE input files.

but for run_process_visc_layer.sh, the output files' layer id is counted from surface (1) and increase downwardly, and it can only process layers belonging to the top cpus. (This is a legacy code.)
