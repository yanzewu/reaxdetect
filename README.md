# ReaxDetect - Reaction Analysis of Atomistic Simulation Trajectory

ReaxDetect aims at extracting time evolution of reaction and species information from atomistic simulation (includes ReaxFF, Ab-inito MD, Quantum Chemistry, etc.). ReaxDetect detects molecule fragment based on connectivity and reactions based on connectivity-change. The species are outputed in Canonicalized ``SMILES`` format. By default, ReaxDetect accepts ReaxFF trajectory as input, but an alternative trajectory parser can be easily written and integrated.

## Example

Analysis for a system with volume=512000 A^3 and frame interval=1.0 ps:

    $ reaxdetect -v 512000 -t 1.0 trajectory.trj

Alternatively:
    
    $ echo '{"volume":512000, "t":1.0}' > config.json
    $ reaxdetect -c config.json trajectory.trj

## Installation

Prerequisties: C++14 supported compiler (*nix); Visual Studio 2017 (Windows).

Install command (*nix):

    make

Windows:

    msbuild /p:Configuration=Release

## Command-line Options:
SYNOPSIS:

	reaxdetect [OPTIONS...] FILE

OPTIONS

- -c, --config [CONFIG]: Load external configuration file. If not specified, the name is 'reacdetect.ini'.
- -b, --RecognizeBegin [BEGIN]: Integer, begin frame number of trajectory reading. Default is 0.
- -s, --RecognizeInterval [BEGIN]: Integer, step of trajectory reading.
- -m, --RecognizeLimit [BEGIN]: Integer, maximum frame number of trajectory reading, Default is -1 (end).
- -t, --timestep [TIMESTEP]: Float, timestep of trajectory file. Default is 0.0 and reading from trajectory.
- -v, --volume [VOLUME]: Volume of system. Default is 1.0.
- -h, --help: Display help.
--version: Display version.

FILE

By default, trajectory file is the control file output from ``LAMMPS`` command:

    pair_style reax/c [control_file]

In the control_file, ``atom_info`` and ``bond_info`` must be set to 1, ``atom_forces``, ``atom_velocities`` and ``angle_info`` must be set to 0.

## Extended Configuration

Configuration file is in json format, with options below:

BondOrderCutoff: Filename of cutoff description file. Use "default" to disable. Format of cutoff description file is:
    
    A-B c0,c1,c2,...

Where A, B is atom name, and A must be __ahead__ of B in periodic table. c0, c1, c2 are step cutoff for bond order 0, 1, 2,... sequentially. Note: 0 is a valid value for bond order, representing weak connections.

If a bond is not in cutoff description, default value is used.

BondOrderCutoffDefault: Cutoff for each type of bond order if not specified in ``BondOrderCutoff``. Format is same (c0,c1,c2,...).

CountBondOrder:     Frame step for bond order statistics. Set '0' to disable. The bond orders will be written in '[name]_bondorder.csv'.

FrameBufferSize:    Same as '-b' option.

ReadAtomPos:        0/1, whether parsing atom positions or simply skipping them in ReaxFF trajectory.


## Output

All files are in csv format.

- Dump file (\*\_full\_dump.csv, \*\_full\_reac.csv): Species and reaction frequency in each frame.
- Report file (\*\_full_report.csv):     Basic analysis of trajectory read, including species lifetimes.
- Bondorder file (\*\_boc.csv):          Bond order of each type of bond. Based on row. Only output when ``CountBondOrder`` != 0.

## Alternative Trajectory Reader

You may write your own trajectory reader based on the interface of ``TrajReader`` at [trajectory.h](reaxdetect/trajectory.h). Details are described in the file.

## Citing

If you use (part of) this program for academic research, please consider citing [].