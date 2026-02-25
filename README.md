# 2025 Testbeam Analysis

This repository contains all code and utilities developed for the 2025 ENUBET testbeam data analysis, in particolar to share code to obtain the timing resolution for our detectors. 

---

## Contents
- C++ scripts for timing and amplitude extraction tools
- Fitting functions for Polya and Gaussian distributions
- ...?

---

## Collaboration

This repository is designed for active collaboration. Feel free to:
- Push experimental branches and test new approaches
- Add notes or comments explaining new methods or assumptions (I'll try to put them at the beginning of each script)
- Suggest improvements to existing scripts

---

## Requirements
- ROOT (latest stable version recommended)
- C++17 or newer

---

## Notes

This repository is evolving continuously as new insights emerge from the testbeam data. Each user is encouraged to document changes and observations directly in the code or via commit messages to maintain clarity and traceability!


---
## Event Building(november test beam)


- 1st part: read the trigger file with bin2root_trigger.cpp this file read the binary file and converts it in a ROOT file.
- 2nd part: read the binary file coming from the SAMPIC with s2root_time_corr.cpp here we have the conversion from bin to root file
- 3rd part: merge of trigger and data using eventbuilder_november.cpp here there is the coupling of the MCP data with the SRS
- (4th part) before the final part: time ordering of the picose event with timeorder.cpp
- last part: evbuilder_november_picosec.cpp reconstruction completed!!!
- NEXT JOB: GEM data for the tracks to add using the anaRun.root file
  
