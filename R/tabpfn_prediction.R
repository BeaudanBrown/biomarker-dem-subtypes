# Loaded inside the pinned GPU R subprocess only. Sourcing defines functions;
# it does not initialize Python, acquire weights, or perform inference.
tabpfn_predict_fold <- function(x, y, newx, config, seed) {
  checkpoint <- normalizePath(config$checkpoint, mustWork = TRUE)
  if (digest::digest(file = checkpoint, algo = "sha256") != config$sha256) {
    stop("TabPFN checkpoint checksum mismatch.")
  }
  Sys.setenv(
    OMP_NUM_THREADS = config$threads,
    OPENBLAS_NUM_THREADS = 1,
    MKL_NUM_THREADS = 1,
    HF_HUB_OFFLINE = "1",
    HF_HUB_DISABLE_TELEMETRY = "1"
  )
  tabpfn::tabpfn_initialize()
  if (
    as.character(utils::packageVersion("tabpfn")) != "0.4.0" ||
      reticulate::py_config()$python != Sys.getenv("RETICULATE_PYTHON")
  ) {
    stop("Unexpected R TabPFN version or Python interpreter.")
  }
  metadata <- reticulate::import("importlib.metadata")
  torch <- reticulate::import("torch")
  if (
    metadata$version("tabpfn") != "9.1.0" ||
      is.null(torch$version$hip) ||
      !torch$cuda$is_available() ||
      config$device != "cuda:0"
  ) {
    stop("Expected local HIP GPU runtime; no CPU/hosted fallback.")
  }
  properties <- torch$cuda$get_device_properties(0L)
  if (
    strsplit(properties$gcnArchName, ":", fixed = TRUE)[[1]][1] != "gfx1030"
  ) {
    stop("Expected the verified gfx1030 GPU.")
  }
  torch$set_num_threads(as.integer(config$threads))
  torch$cuda$empty_cache()
  torch$cuda$reset_peak_memory_stats()
  control <- tabpfn::control_tab_pfn(
    model_path = checkpoint,
    device = config$device,
    random_state = as.integer(seed),
    n_preprocessing_jobs = 1L,
    ignore_pretraining_limits = FALSE,
    inference_precision = "auto",
    fit_mode = "fit_preprocessors",
    memory_saving_mode = "auto"
  )
  model <- tabpfn::tab_pfn(
    x,
    factor(y, levels = 0:1),
    num_estimators = config$members,
    softmax_temperature = 0.9,
    average_before_softmax = FALSE,
    balance_probabilities = FALSE,
    version = NULL,
    control = control
  )
  predictions <- predict(model, newx, type = "prob")
  if (!".pred_1" %in% names(predictions))
    stop("Missing positive-class probability .pred_1.")
  predictions <- as.numeric(predictions[[".pred_1"]])
  if (
    length(predictions) != nrow(newx) ||
      any(!is.finite(predictions)) ||
      any(predictions < 0 | predictions > 1)
  )
    stop("Invalid TabPFN probabilities.")
  torch$cuda$synchronize()
  fitted <- model$fit
  if (
    fitted$n_estimators_ != config$members ||
      fitted$n_train_samples_ != nrow(x) ||
      !is.null(fitted$forced_inference_dtype_) ||
      !fitted$use_autocast_
  ) {
    stop("Unexpected effective TabPFN configuration.")
  }
  builtins <- reticulate::import_builtins()
  devices <- vapply(
    fitted$devices_,
    function(device) builtins$str(device),
    character(1)
  )
  if (!identical(unname(devices), config$device))
    stop("Unexpected fitted device.")
  allocated <- torch$cuda$max_memory_allocated() / 2^20
  if (allocated <= 0) stop("No GPU allocations recorded.")
  diagnostics <- list(
    device = devices,
    gpu = properties$name,
    members = fitted$n_estimators_,
    training_rows = fitted$n_train_samples_,
    autocast = fitted$use_autocast_,
    parameters = builtins$repr(fitted$get_params(deep = FALSE)),
    member_configs = builtins$repr(fitted$ensemble_configs_),
    peak_allocated_mib = allocated,
    peak_reserved_mib = torch$cuda$max_memory_reserved() / 2^20,
    python = reticulate::py_config()$python,
    r_version = as.character(getRversion()),
    r_tabpfn = as.character(utils::packageVersion("tabpfn")),
    python_tabpfn = metadata$version("tabpfn"),
    torch = metadata$version("torch"),
    hip = torch$version$hip,
    checkpoint_sha256 = config$sha256
  )
  # Return only serializable R values; child exit releases the GPU model/cache.
  list(predictions = predictions, diagnostics = diagnostics)
}
