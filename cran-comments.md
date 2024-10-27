## Test environments
* local Manjaro Linux 24.1.1 install, R 4.4.1
* win-builder (release, oldrelease, devel)
* mac-builder
* rhub (platforms: linux, macos, macos-arm64, windows)


## R CMD check results

### local
0 errors | 0 warnings | 0 notes

### R-hub
0 errors | 0 warnings | 1 note

* checking Rd cross-references ... NOTE
Packages unavailable to check Rd xrefs: ‘DiscreteFDR’, ‘FDX’

- rhub's fault; these packages are obviously not installed on their
  VMs/containers

### win-builder
0 errors | 0 warnings | 0 notes

### mac-builder
0 errors | 0 warnings | 0 notes
