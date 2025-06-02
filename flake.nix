{
  description = "A basic flake with a shell";
  outputs = { self, nixpkgs, flake-utils, pre-commit-hooks, ... }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
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

        devShells.default = pkgs.mkShell {
          inherit (self.checks.${system}.pre-commit-check) shellHook;
          env.R_LIBS_USER="./.Rlib";
          buildInputs = [
            pkgs.bashInteractive
            self.checks.${system}.pre-commit-check.enabledPackages
          ];
          packages = with pkgs;
            [
              R
              quarto
              air-formatter
            ] ++ (with rPackages; [
              languageserver
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
              nnet
              cvAUC
              gbm
              e1071
              class
              patchwork
              targets
              usethis
              visNetwork
              qs2
              styler
              tarchetypes
              bartMachine
            ]);
        };
      }
    );

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
    pre-commit-hooks.url = "github:cachix/git-hooks.nix";
    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.systems.follows = "systems";
    };
  };
}
