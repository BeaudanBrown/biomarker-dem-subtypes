# User-run GPU subprocess entrypoint, never part of ordinary checks.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("Expected input, output and prediction-adapter paths.")
source(args[3], local = TRUE)
input <- readRDS(args[1])
result <- do.call(tabpfn_predict_fold, input)
saveRDS(result, args[2])
Sys.chmod(args[2], "0600")
