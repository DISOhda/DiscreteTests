## Test environments
* local Windows 10 Pro 23H2 install, R 4.3.3
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
0 errors | 0 warnings | 1 note
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Florian Junge <diso.fbmn@h-da.de>'

New submission

Possibly misspelled words in DESCRIPTION:
  vectorised (11:23)
  
- ever heard of BRITISH english? (check seems to ignore that "Lang" field in
  DESCRIPTION is set to "en-GB")
