This a critical update that fixes a bug in `mann_whitney_test_pv()` that
consistently produces wrong p-values if normal approximation is used.

## Test environments
* local Manjaro Linux 26.1.1 install, Kernel 7.1.9, R 4.6.1
* win-builder (release, oldrelease, devel)
* mac-builder (release, devel)
* rhub (platforms: linux, m1-san, macos, macos-arm64, windows, valgrind), see
  https://github.com/DISOhda/DiscreteTests/actions/runs/33079053199


## R CMD check results

### local
0 errors | 0 warnings | 0 notes

### R-hub
0 errors | 0 warnings | 0 notes

### win-builder
0 errors | 0 warnings | 0 notes

### mac-builder
0 errors | 0 warnings | 0 notes
