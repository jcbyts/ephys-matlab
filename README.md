# ephys-matlab
Matlab tools for electrophysiology experiments collected with [PLDAPS](https://github.com/huklab/PLDAPS/tree/openreception) and [Open-Ephys](https://github.com/open-ephys/plugin-GUI)

### General idea:
ephys-matlab provides a set of tools for essential components of all visual neuroscience experiments:
 1. importing 
 	* synchronizing visual stimuli, behavior, and electrophysiology
 	* convert from different file-types
 	* handling channel maps from multiple recording arrays/headstages
 2. pre-processing
 	* local field potential (LFP), filtering, line-noise removal, phase correction
 	* saccadic eye movements detection
 	* spike sorting via KiloSort or waveform clustering
3. reverse correlation
	* build design matrix

caveats: This is a work in progress and will gain new features as they become necessary. This package will not be supported for public consumption.

PLDAPS manages stimulus presentation and behavioral monitoring. Open-Ephys 

## Importing a session

Make sure that all the paths are added
```matlab
addEphysMatlab
```