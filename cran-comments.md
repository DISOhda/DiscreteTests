## Test environments
* local Manjaro Linux 26.0.2 install, R 4.5.2
* win-builder (release, oldrelease, devel)
* mac-builder
* rhub (platforms: linux, m1-san, macos-arm64, windows, lto, valgrind), see
  https://github.com/DISOhda/DiscreteTests/actions/runs/21975650989


## R CMD check results

### local
0 errors | 0 warnings | 0 notes

### R-hub
0 errors | 0 warnings | 1 note

* checking Rd cross-references ... NOTE
Packages unavailable to check Rd xrefs: ‘DiscreteFDR’, ‘FDX’

  - rhub's fault; these packages are obviously not installed on their
    VMs/containers
  - does only occur on some configurations

### win-builder
0 errors | 0 warnings | 1 note

* checking DESCRIPTION meta-information ... NOTE
  Author field differs from that derived from Authors@R
  
  - must be false alarm, there is no additional Authors field in DESCRIPTION
    file
  - happens only for oldrelease checks, not for release or devel

### mac-builder
0 errors | 0 warnings | 0 notes
