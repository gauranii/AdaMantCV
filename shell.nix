# R with the three packages the scripts need, pinned to one nixpkgs commit.
let pkgs = import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/f3aa7d3d72c2bc66e0f30f4262d7ec2e33c9eab2.tar.gz") {};
in pkgs.mkShell { buildInputs = [ (pkgs.rWrapper.override { packages = with pkgs.rPackages; [ ggplot2 MASS glmnet ]; }) ]; }
