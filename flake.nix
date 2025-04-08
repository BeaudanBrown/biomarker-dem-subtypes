{
  description = "A basic flake with a shell";
  outputs = { nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
        {
        devShells.default = pkgs.mkShell {
          env.R_LIBS_USER="./.Rlib";
          buildInputs = [ pkgs.bashInteractive ];
          packages = with pkgs;
            [
              R
              quarto
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
            ]);
        };
      }
    );

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    systems.url = "github:nix-systems/default";
    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.systems.follows = "systems";
    };
  };
}
