# ReaxDetect Manual 

## Installation

Requirement: C++14 supported compiler; GCC 5.0+ recommended.

Install command:

    make

## Usage: 

	reaxdetect [OPTIONS...] FILE

### OPTIONS

-b [buffer_size]:   How many frames used to scan reaction. Default is 2. If large than 2, will eliminate reversed reactions after 'buffer_size-1' frames after a recognition.

-c [config_file]:   Load external configuration file. If not specified, the name is 'reacdetect.ini'.

--config:           Write default configuration to file 'reacdetect.ini'.

--dump=[dumpoption]:  Control raw output (reactions and species in each frame). dumpoption can be:
- nodump: Do not dump;
- full: Dump all.

-f [step]:          How many frames as a step in reading trajectory file.

-h, --help:         Display help.

-l [limit]:         Maximum number of frame recognized. 

-s [sampeinterval]: Number of frames used in sample.

-t [timestep]:      Timestep of trajectory file. Default is reading from trajectory.

-v [volume]:        Volume of system. Default is 1.0.

--version:          Display version.

### FILE

The trajectory file is specified from lammps command:

    pair_style reax/c [control_file]

In the control_file, "atom_info" and "bond_info" must be set to 1. While "atom_forces", "atom_velocities" and "angle_info" must be set to 0.


### Configuration file

BondOrderCutoff:    Filename of cutoff description file. Use "default" to use default values. Format of cutoff description file:
    
    A-B c0,c1,c2,...

Where A, B is atom name, and A must be ahead of B in periodic table. c0, c1, c2 are cutoff for different bond order.
If a bond is not in cutoff description, default value is used.

BondOrderCutoffDefault: Cutoff for each type of bond order if not specified in BondOrderCutoff. Format is same as above.

CountBondOrder:     Frame step for bond order statistics. Set '0' to disable. The bond orders will be written in '[name]_bondorder.csv'.

FrameBufferSize:    Same as '-b' option.

ReadAtomPos:        true/false. For future use.

RecognizeBegin:     Begin frame number of trajectory analyzation.

RecognizeInterval:  Same as '-f' option.

RecognizeLimit:     Same as '-l' option.

SampleMethod:       How to get samples in trajectory. Currently only "fixint" is available.

SampleInterval:     Same as '-s' option.

SampleRange:        Range for sampling. If different from 'SampleInterval', take only 'SampleRange' frames to average in sampling.


## Output

All files are in csv format.

Raw frequency file (\*\_rawfreq.csv):  Data for kinetic calculation based on power law.

Sample file (\*\_sample.csv):          Sample of reaction trajectory (how many species and reactions).

Report file (\*\_full_report.csv):     Basic data in analyzation. Life of species.

Dump file (\*\_full\_dump.csv, \*\_full\_reac.csv): Species and reaction number in each frame, respectively. Only output when dump=full.

Bondorder file (\*\_boc.csv):          Bond order of each type of bond. Based on row. Only 


## Example

Analyzation for a system with volume=512000 A^3 and frame interval=1.0 ps:

    $ reaxdetect --config
    $ reaxdetect -v 512000 -t 1.0 trajectory.trj