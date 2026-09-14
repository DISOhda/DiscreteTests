This a critical update that fixes a bug in `wilcox_test_pv()` that
can produces wrong p-values if normal approximation on samples with many
zero-differences is used. Sometimes, NAs are produced because of incorrectly
computed variances (could become negative!)

## Test environments
* local Manjaro Linux 26.1.2 install, Kernel 7.1.13, R 4.6.1
* win-builder (release, oldrelease, devel)
* no mac-builder tests, since it was unavailable; substituted with rhub tests
* rhub (platforms: linux, m1-san, macos, macos-arm64, windows, valgrind), see
  https://github.com/DISOhda/DiscreteTests/actions/runs/34824959588/


## R CMD check results

### local
0 errors | 0 warnings | 0 notes

### win-builder
0 errors | 0 warnings | 0 notes

### R-hub
0 errors | 0 warnings | 0 notes
