# ReaxDetect Manual 

## Installation

Requires c++14 support. Command:

    make 

## Usage: 

	reaxdetect [-b buffersize] [-c config] [--config] [--dump=nodump/full] [-f step] [-h,--help] [-l limit] [-s samplerange] [-t timestep] [-v volume] [--version] trajectory

### Argument description:

-b buffer_size: How many frames used to scan reaction. Must be >= 2. Default is 2.

-c config:      Load external configuration file.

--config:       Write default configuration file.

--dump=dumpoption:  Control raw output (reactions and species in each frame). nodump: Do not dump; full: Dump all.

-f step:        Frame skipped.

-h, --help:     Display help.

-l limit:       Maximum frame recognized.

-s samplerange: Number of frames used in sample.

-t timestep:    Timestep of trajectory file. Default is reading from trajectory.

-v volume:      Volume of system. Default is 1.

--version:      Display version.


### Configuration file

BondOrderCutoff:    Filename of cutoff description file. Use "default" to use default values. Format of cutoff description file:
    
    A-B c0,c1,c2,...

Where A, B is atom name, and A must be ahead of B in periodic table. c0, c1, c2 are cutoff for different bond order.
If a bond is not in cutoff description, default value is used.

BondOrderCutoffDefault: Cutoff for each type of bond order if not specified in BondOrderCutoff. Format is same as above.

CountBondOrder:     true/false. Specify bond order histogram output.

FrameBufferSize:    How many frames used to scan reaction. Must be >= 2.

ReadAtomPos:        Read atom position in trajectory.

RecognizeBegin:     Begin frame number of trajectory analyzation.

RecognizeInterval:  Step of trajectory analyzation.

RecognizeLimit:     Maximum frame number of trajectory analyzation.

SampleMethod:       How to get samples in trajectory. "fixint": Fixed interval.

SampleInterval:     Interval for sampling.

SampleRange:        Range for sampling.


## Output

All files are in csv format.

Raw frequency file (_rawfreq.csv):  Data for kinetic calculation based on power law.

Sample file (_sample.csv):          Sample of reaction trajectory (how many species and reactions).

Report file (_full_report.csv):     Basic data in analyzation. Life of species.

Dump file (_full_dump.csv, _full_reac.csv): Species and reaction number in each frame, respectively. Only output when dump=full.

Bondorder file (_boc.csv):          Bond order of each type of bond. Based on row.