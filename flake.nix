{
  description = "AdaMantCV: code for Gauran, Xue, Chen, Ombao, Yu (2022), Front. Neurosci. 16:836100";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/f3aa7d3d72c2bc66e0f30f4262d7ec2e33c9eab2";

  outputs =
    { self, nixpkgs }:
    let
      lib = nixpkgs.lib;
      systems = [
        "aarch64-darwin"
        "x86_64-darwin"
        "aarch64-linux"
        "x86_64-linux"
      ];
      forEachSystem = f: lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});

      # R with the packages the scripts load.
      rEnv =
        pkgs:
        pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            ggplot2
            MASS
            glmnet
          ];
        };

      # Python with the packages of python/adamant; numba compiles the random
      # number generator and the loops that follow R's order of operations.
      pyEnv =
        pkgs:
        pkgs.python3.withPackages (
          ps: with ps; [
            numpy
            scipy
            pandas
            matplotlib
            numba
          ]
        );

      # The BLAS and LAPACK library R is linked to; python/adamant/rbase.py
      # calls it for its matrix products and decompositions.
      blasLibrary =
        pkgs: "${pkgs.openblas}/lib/libopenblas${pkgs.stdenv.hostPlatform.extensions.sharedLibrary}";

      # One command per table and figure of the paper.
      scripts = {
        table1 = "table1_type1_error";
        table2 = "table2_power_fixed";
        table3 = "table3_power_random";
        table4 = "table4_eeg_pvalues";
        fig1 = "fig1_liver_cv";
        fig2 = "fig2_liver_heat";
        fig3 = "fig3_lambda_hist";
        fig5 = "fig5_eeg_heat";
        fig6 = "fig6_eeg_heat";
        fig7 = "fig7_power_gcv_decreasing";
        fig8 = "fig8_power_gcv_increasing";
      };

      # Run the scripts from the flake source; results/ and figures/ are
      # written to the directory the command is started in. The scripts run
      # their replications in parallel R processes (--cores), so each process
      # gets one BLAS thread.
      runner =
        pkgs: name: names:
        pkgs.writeShellScript "adamantcv-${name}" ''
          set -euo pipefail
          export ADAMANT_OUTPUT="''${ADAMANT_OUTPUT:-$PWD}"
          export OPENBLAS_NUM_THREADS=1
          cd ${self}
          for script in ${lib.concatStringsSep " " names}; do
            ${rEnv pkgs}/bin/Rscript "scripts/$script.R" "$@"
          done
        '';
      app = pkgs: name: names: {
        type = "app";
        program = "${runner pkgs name names}";
        meta.description = "Build ${name} of the paper";
      };

      # The same commands with the Python scripts, as py-<name>.
      pyRunner =
        pkgs: name: names:
        pkgs.writeShellScript "adamantcv-py-${name}" ''
          set -euo pipefail
          export ADAMANT_OUTPUT="''${ADAMANT_OUTPUT:-$PWD}"
          export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
          export ADAMANT_BLAS=${blasLibrary pkgs}
          cd ${self}
          for script in ${lib.concatStringsSep " " names}; do
            ${pyEnv pkgs}/bin/python "python/scripts/$script.py" "$@"
          done
        '';
      pyApp = pkgs: name: names: {
        type = "app";
        program = "${pyRunner pkgs name names}";
        meta.description = "Build ${name} of the paper with the Python scripts";
      };
    in
    {
      devShells = forEachSystem (pkgs: {
        default = pkgs.mkShell {
          packages = [
            (rEnv pkgs)
            (pyEnv pkgs)
          ];
          OPENBLAS_NUM_THREADS = 1;
          VECLIB_MAXIMUM_THREADS = 1;
          ADAMANT_BLAS = blasLibrary pkgs;
        };
      });

      apps = forEachSystem (
        pkgs:
        lib.mapAttrs (name: script: app pkgs name [ script ]) scripts
        // lib.mapAttrs' (name: script: lib.nameValuePair "py-${name}" (pyApp pkgs name [ script ])) scripts
        // {
          all = app pkgs "all" (lib.attrValues scripts);
          py-all = pyApp pkgs "all" (lib.attrValues scripts);
        }
      );

      # Build every table and figure with the default settings and compare
      # each with the paper; the result is the comparison (checks/RESULTS.md).
      checks = forEachSystem (pkgs: {
        outputs =
          pkgs.runCommand "adamantcv-outputs"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                pkgs.python3
              ];
            }
            ''
              export HOME=$TMPDIR ADAMANT_OUTPUT=$TMPDIR
              cd ${self}
              bash checks/run_checks.sh $TMPDIR/RESULTS.md || { cat $TMPDIR/RESULTS.md; exit 1; }
              cp $TMPDIR/RESULTS.md $out
            '';

        # Build every table and figure with the R and with the Python
        # scripts and compare the two, and compare the random numbers,
        # pgamma and the cross-validated ridge regression of the Python
        # code with R; the result is the comparison (checks/PYTHON.md, with
        # the stored comparison of the reruns).
        python =
          pkgs.runCommand "adamantcv-python"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                (pyEnv pkgs)
              ];
              ADAMANT_BLAS = blasLibrary pkgs;
            }
            ''
              export HOME=$TMPDIR
              cd ${self}
              bash checks/check_python.sh archive $TMPDIR/work $TMPDIR/PYTHON.md $TMPDIR \
                || { cat $TMPDIR/PYTHON.md; exit 1; }
              python python/tests/test_rng.py
              python python/tests/test_rmath.py
              python python/tests/test_glmnet.py quick
              cp $TMPDIR/PYTHON.md $out
            '';
      });

      formatter = forEachSystem (pkgs: pkgs.nixfmt-rfc-style);
    };
}
