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
    in
    {
      devShells = forEachSystem (pkgs: {
        default = pkgs.mkShell {
          packages = [
            (rEnv pkgs)
            pkgs.python3
          ];
          OPENBLAS_NUM_THREADS = 1;
        };
      });

      apps = forEachSystem (
        pkgs:
        lib.mapAttrs (name: script: app pkgs name [ script ]) scripts
        // {
          all = app pkgs "all" (lib.attrValues scripts);
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
      });

      formatter = forEachSystem (pkgs: pkgs.nixfmt-rfc-style);
    };
}
