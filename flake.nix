{
  description = "A basic flake with a shell";
  outputs =
    {
      self,
      nixpkgs,
      nixpkgs-tabpfn,
      flake-utils,
      pre-commit-hooks,
      pyproject-nix,
      uv2nix,
      pyproject-build-systems,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        # Independent pin: do not upgrade the existing analysis R/packages.
        gpuPkgs = nixpkgs-tabpfn.legacyPackages.${system};
        tabpfnGpu = import ./environments/tabpfn-gpu {
          pkgs = gpuPkgs;
          inherit pyproject-nix uv2nix pyproject-build-systems;
        };
        tabpfnR = import ./environments/tabpfn-gpu/r.nix {
          pkgs = gpuPkgs;
          pythonEnvironment = tabpfnGpu.environment;
        };
      in
      {
        checks = {
          pre-commit-check = pre-commit-hooks.lib.${system}.run {
            src = ./.;
            hooks = {
              air-fmt = {
                enable = true;
                entry = "air format";
                files = ".*\.[rR]$";
              };
            };
          };
        };

        devShells = {
          default = pkgs.mkShell {
            inherit (self.checks.${system}.pre-commit-check) shellHook;
            env.R_LIBS_USER = "./.Rlib";
            buildInputs = [
              pkgs.bashInteractive
              self.checks.${system}.pre-commit-check.enabledPackages
            ];
            packages =
              with pkgs;
              [
                R
                quarto
                air-formatter
              ]
              ++ (with rPackages; [
                languageserver
                processx
                withr
                dotenv
                readxl
                performance
                lme4
                data_table
                tidyverse
                SuperLearner
                gtsummary
                ggthemes
                vimp
                xgboost
                glmnet
                ranger
                randomForest
                earth
                gam
                pROC
                origami
                future
                future_apply
                logistf
                rms
                psych
                bayesplot
                knitr
                cardx
                corrr
                quarto
                arm
                cowplot
                plotly
                htmlwidgets
                htmltools
                nnet
                cvAUC
                gbm
                e1071
                class
                patchwork
                targets
                crew
                usethis
                visNetwork
                qs2
                styler
                tarchetypes
                dbarts
                tmle
                kableExtra
              ]);
          };
        }
        // pkgs.lib.optionalAttrs (system == "x86_64-linux") {
          tabpfn-gpu = self.devShells.${system}.default.overrideAttrs (_old: {
            ADDF_TABPFN_RSCRIPT = "${tabpfnR}/bin/Rscript";
          });
        };
      }
    );

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    nixpkgs-tabpfn.url = "github:NixOS/nixpkgs/nixos-26.05";
    pyproject-nix = {
      url = "github:pyproject-nix/pyproject.nix";
      inputs.nixpkgs.follows = "nixpkgs-tabpfn";
    };
    uv2nix = {
      url = "github:pyproject-nix/uv2nix";
      inputs.nixpkgs.follows = "nixpkgs-tabpfn";
      inputs.pyproject-nix.follows = "pyproject-nix";
    };
    pyproject-build-systems = {
      url = "github:pyproject-nix/build-system-pkgs";
      inputs.nixpkgs.follows = "nixpkgs-tabpfn";
      inputs.pyproject-nix.follows = "pyproject-nix";
      inputs.uv2nix.follows = "uv2nix";
    };
    systems.url = "github:nix-systems/default";
    pre-commit-hooks.url = "github:cachix/git-hooks.nix";
    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.systems.follows = "systems";
    };
  };
}
