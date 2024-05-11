## Test environments
* local Windows 10 Pro 23H2 install, R 4.4.0
* win-builder (release, oldrelease, devel)
* R-hub (platforms: linux, macos, macos-arm64, windows)


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
