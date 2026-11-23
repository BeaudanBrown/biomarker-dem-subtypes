{
  pkgs,
  pyproject-nix,
  uv2nix,
  pyproject-build-systems,
}:
let
  # patchelf 0.15.2 corrupts relocated PT_PHDR in AMD's libhipblaslt.
  # The newer patchelf already pinned by nixpkgs preserves its mapped headers.
  # Scope this tool change to GPU packaging; leave CPU/R packages untouched.
  workspace = uv2nix.lib.workspace.loadWorkspace { workspaceRoot = ./.; };
  overlay = workspace.mkPyprojectOverlay { sourcePreference = "wheel"; };
  pythonSet =
    (pkgs.callPackage pyproject-nix.build.packages {
      python = pkgs.python312;
    }).overrideScope
      (
        pkgs.lib.composeManyExtensions [
          pyproject-build-systems.overlays.wheel
          overlay
          (
            _final: prev:
            pkgs.lib.genAttrs
              [
                "rocm-sdk-core"
                "rocm-sdk-libraries"
                "rocm-sdk-device-gfx1030"
                "amd-torch-device-gfx1030"
                "torch"
                "triton"
              ]
              (
                name:
                prev.${name}.overrideAttrs (old: {
                  nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ [ pkgs.patchelfUnstable ];
                  preFixup = ''
                    export PATH="${pkgs.patchelfUnstable}/bin:$PATH"
                  ''
                  + (old.preFixup or "");
                })
              )
          )
          (final: prev: {
            rocm-sdk-core = prev.rocm-sdk-core.overrideAttrs (old: {
              buildInputs = (old.buildInputs or [ ]) ++ [ pkgs.ocl-icd ];
              # The SDK bundles debugger executables for other Python ABIs.
              # This environment supports only 3.12; don't pull four interpreters
              # and legacy libcrypt into the runtime for unused debugger variants.
              postInstall = (old.postInstall or "") + ''
                rm -f "$out/${final.python.sitePackages}/_rocm_sdk_core/bin/rocgdb-py3."{10,11,13,14}
              '';
            });
            rocm-sdk-libraries = prev.rocm-sdk-libraries.overrideAttrs (old: {
              buildInputs = (old.buildInputs or [ ]) ++ [ final.rocm-sdk-core ];
              # pip merges the device payload into the same Python namespace.
              # dladdr/RPATH resolve native libraries to their own store output,
              # so those libraries also need their adjacent kernel/data files.
              postInstall = (old.postInstall or "") + ''
                cp -as --no-clobber ${final.rocm-sdk-device-gfx1030}/${final.python.sitePackages}/_rocm_sdk_libraries/. "$out/${final.python.sitePackages}/_rocm_sdk_libraries/"
              '';
              preFixup = (old.preFixup or "") + ''
                addAutoPatchelfSearchPath ${final.rocm-sdk-core}/${final.python.sitePackages}/_rocm_sdk_core/lib
              '';
            });
            torch = prev.torch.overrideAttrs (old: {
              postInstall = (old.postInstall or "") + ''
                mkdir -p "$out/${final.python.sitePackages}/torch/.kpack"
                ln -s ${final.amd-torch-device-gfx1030}/${final.python.sitePackages}/torch/.kpack/torch_gfx1030.kpack "$out/${final.python.sitePackages}/torch/.kpack/torch_gfx1030.kpack"
              '';
              buildInputs = (old.buildInputs or [ ]) ++ [
                final.rocm-sdk-core
                final.rocm-sdk-libraries
                pkgs.xz
              ];
              preFixup = (old.preFixup or "") + ''
                addAutoPatchelfSearchPath ${final.rocm-sdk-core}/${final.python.sitePackages}/_rocm_sdk_core/lib
                addAutoPatchelfSearchPath ${final.rocm-sdk-libraries}/${final.python.sitePackages}/_rocm_sdk_libraries/lib
              '';
            });
            triton = prev.triton.overrideAttrs (old: {
              buildInputs = (old.buildInputs or [ ]) ++ [
                final.rocm-sdk-core
                final.rocm-sdk-libraries
                pkgs.zlib
              ];
              preFixup = (old.preFixup or "") + ''
                addAutoPatchelfSearchPath ${final.rocm-sdk-core}/${final.python.sitePackages}/_rocm_sdk_core/lib
                addAutoPatchelfSearchPath ${final.rocm-sdk-libraries}/${final.python.sitePackages}/_rocm_sdk_libraries/lib
              '';
            });
            rocm = prev.rocm.overrideAttrs (old: {
              nativeBuildInputs =
                (old.nativeBuildInputs or [ ])
                ++ final.resolveBuildSystem {
                  setuptools = [ ];
                };
            });
          })
        ]
      );
  environment = pythonSet.mkVirtualEnv "tabpfn-gfx1030-environment" workspace.deps.default;
in
{
  inherit environment;
}
