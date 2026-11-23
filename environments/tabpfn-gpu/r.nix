# Adapted from ~/collab/tabpfn-test/environments/tabpfn/r.nix.
# Keep this R/Python stack separate from the existing ADDF analysis runtime.
{ pkgs, pythonEnvironment }:
let
  tabpfn = pkgs.rPackages.buildRPackage {
    name = "tabpfn";
    version = "0.4.0";
    src = pkgs.fetchurl {
      urls = [
        "https://cran.r-project.org/src/contrib/tabpfn_0.4.0.tar.gz"
        "https://cran.r-project.org/src/contrib/Archive/tabpfn/tabpfn_0.4.0.tar.gz"
      ];
      hash = "sha256-4/soDQdKNPGcWuiyUJtqBYXS70tk6r6OXApOzY/wD00=";
    };
    propagatedBuildInputs = with pkgs.rPackages; [
      cli
      dplyr
      generics
      hardhat
      jsonlite
      purrr
      reticulate
      rlang
      tibble
    ];
    # Upstream tests may provision weights/Python; never run them at build time.
    doCheck = false;
  };
  baseR = pkgs.rWrapper.override {
    packages = [
      tabpfn
      pkgs.rPackages.digest
    ];
  };
in
pkgs.symlinkJoin {
  name = "addf-tabpfn-gpu-r";
  paths = [ baseR ];
  nativeBuildInputs = [ pkgs.makeWrapper ];
  postBuild = ''
    for program in R Rscript; do
      wrapProgram "$out/bin/$program" \
        --unset R_HOME --unset R_LIBS --set R_LIBS_USER /nonexistent --unset R_LIBS_SITE \
        --unset PYTHONPATH --unset PYTHONHOME --unset TABPFN_TOKEN \
        --unset LD_LIBRARY_PATH --unset LD_PRELOAD \
        --unset ROCM_PATH --unset HIP_PATH --unset HSA_OVERRIDE_GFX_VERSION \
        --set RETICULATE_PYTHON ${pythonEnvironment}/bin/python \
        --set RETICULATE_USE_MANAGED_VENV no \
        --set PYTHONNOUSERSITE 1 --set UV_PYTHON_DOWNLOADS never \
        --set HF_HUB_OFFLINE 1 --set HF_HUB_DISABLE_TELEMETRY 1
    done
  '';
}
