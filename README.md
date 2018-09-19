# ephys-matlab
Matlab tools for electrophysiology experiments collected with [PLDAPS](https://github.com/huklab/PLDAPS/tree/openreception) and [Open-Ephys](https://github.com/open-ephys/plugin-GUI)

### General idea:
ephys-matlab provides a set of tools for essential components of visual neuroscience experiments:
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

## Meta Data
A central component of ephys-matlab is the use of meta data for tracking what experiments / analyses have been run. In many ways, we use this as if we had a database, but it is not one!! It is just an excel file that keeps track of the

### Automagically tracked features
* Date
* Subject
* Directory
* Time
* Tag
* oe2dat
* StimulusProtocols
* SpikeSorting
* LFpPhaseCorrection

### Manually tracked features
* Weight
* Chamber
* Rig
* Electrode
* Lens
* FlipXEye
* FlipYEye


## Installation
The following git repositories are required for ephys-matlab to import a session.

Repos associated with PLDAPS/Eyelink/OpenEphys:

[PLDAPS](https://github.com/jcbyts/PLDAPS/tree/noparams), 
[PEP](https://github.com/jcbyts/pds-stimuli), 
[edfmex](https://github.com/HukLab/edfmex),
[analysis-tools](https://github.com/open-ephys/analysis-tools)

Repos associated with KiloSort:

[KiloSort](https://github.com/cortex-lab/KiloSort), 
[Spkes](https://github.com/cortex-lab/spikes), 
[SortQuality](https://github.com/cortex-lab/sortingQuality), 
[npy-matlab](https://github.com/kwikteam/npy-matlab)

Once the repositories are on you disk, edit `addEphysMatlab` to appropriately point to each of them and then add them to the matlab path.

`ephys-matlab` uses matlab preferences to manage the location of local data and server data. Set up both directories using the following commands in matlab.

`setpref('EPHYS', 'SERVER_DATA', '/full/path/to/the/data/on/the/server')`

`setpref('EPHYS', 'LOCAL_DATA', '/full/path/to/folder/on/local/disk')`

## Importing a session
The file `import_script.m` provides a template for importing an ephys session collected with PLDAPS / Open-Ephys / Eyelink

In the future it will support other eyetrackers, but for the time being, it expects that every session will contain three types of files
1. Open-Ephys files (.continuous, .spikes, .events) - These are binary files. They contain the raw electrophysiological data and meta data about the session.
2. PLDAPS files (*.PDS) - These are .mat files created by PLDAPS v4 ([noparams](https://github.com/jcbyts/PLDAPS/tree/noparams)). they contain the stimulus and behavioral data, as well as the meta data required to synchronize the Ephys and Display computer clocks
2. Eyelink files (*.edf) - These are a proprietrary filetype from SR research for storying eyelink data.

Make sure that all the paths are added
```matlab
addEphysMatlab
```

### Setup the electrodes / headstages used
"Shank" refers to a particular array or group of electrodes. Ephys-matlab will separate recordings by shank. Shanks will be spike-sorted separately and the processed data will live in separate directories. Use `edit hardware.electrodeFactory` to see the available recording devices. To add a recording device, simply add a new `case` to the switch statement and update it to point to the appropriate files (under `hardware.electrode`, and `hardware.headstage`).

Here is an example of a recording shank entry:

``` matlab
shank{1} = hardware.electrode.Shank2;
shank{1}.headstages{1} = hardware.headstage.intan_RHD2132; % specify what type of headstage was used
shank{1}.name = 'V1'; % name the electrode array
```

`hardware.electrode.(probeName)` are all `probe` objects. You have to create one that matches the probe you are using. The important properties are `channelMap` and `xcoords`, `ycoords`, `zcoords`. Counter-intuitively, `x` and `z` code placement in the chamber and `ycoords` codes for the depth of each electrode site.
The line `shank{1} = hardware.electrode.Shank2;` points to the 32 channel silicone linear array that we have named "shank2".
 ``` matlab
 Shank2 with properties:

            name: 'Shank2'
    manufacturer: 'Atlas'
          design: 'E32-50-S1-L6'
             num: '2015-151'
         xcoords: [32×1 double]
         ycoords: [32×1 double]
         zcoords: [32×1 double]
      channelMap: [1×32 double]
       connector: 'Omnetics36'
        material: 'IrOx'
       impedence: 0.5000
      headstages: {}
```

`hardware.headstage.(headstageName)` refer to the digitizing headstages used with intan/open-ephys. 

For example, we have several 32 channel headstages from intan.
``` matlab
intan_RHD2132 with properties:

            name: 'intan_RHD2132'
    manufacturer: 'intan'
           model: 'RHD2132'
          filter: [1 7500]
    samplingRate: 30000
       connector: 'Omnetics36'
           gains: NaN
      channelMap: [1×36 double]
```

Headstages can apply a channel map to a probe, which will re-order the channels properly. Setting up the headstage and electrode per the manufacturers specs will produce a channel map that accurately orders the channels.

Headstages can also be blank and apply no channel map. This assumes that the channel map specified by the electrode is the correct mapping. This might be useful if using single electrodes, when you already know which channel corresponds to which electrode and there is no particular geometry.

For example: For single electrodes specify a custom channel map using numbers relative to the headstage start. eg., if using ch36 plugged into headstage 2, this should be channel 4 (relative to the start of headstage 2)
```matlab
% list single electrode channels (this can be a vector if > 1 electrode used)
chanMap = 4;
shank{2} = hardware.electrode.customChannelMap(chanMap);
shank{2}.name = 'MtBurrHoleMapping';
```
If you chose hardware.electrode.customChannelMap(chNum), the channel map will be chanMap. That's it. No specifying headstage necessary. Behind the scenes, it is calling `hardware.headstage.blank32`, which applies no channel map.

### run import session
Call `io.importSession()` without any arguments to run a generic import. This will ask you to select a directory where the data from the session live. Then it will prompt for a recording device (shank) and then automatically population the meta table with information from the available data.

## Working with ephys-matlab
If a session has been imported, it will show up as a row in the excel file under `getpref('EPHYS', 'SERVERDATA')`.

For example,

``` matlab
meta = io.getExperimentsAnd('Subject', 'Ellie', 'Chamber', 'V1', 'Lens', 1, 'StimulusProtocols', 'hartleyFF', 'SpikeSorting', 'Kilo');
```



