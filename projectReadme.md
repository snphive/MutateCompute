# Project Title

A collection of scripts and classes to run a family of algorithms (FoldX, Tango, Waltz, Limbo) on protein 
structure files (pdb) and protein sequence files (FASTA).
This includes  

## Getting Started

To run from Pycharm IDE, select Run/Run from toolbar.   
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The current project interpreter: Python 3.7.
The current installed package dependencies are: 
PyYAML 3.13
biopython 1.72
mysql-connector-python 8.0.13
natsort 5.5.0
numpy 1.15.4
parse 1.9.0
parse-type 0.4.2
pip 18.1
protobuf 3.6.1
pydevd 1.4.0
setuptools 40.6.2
six 1.11.0

### Installing

A step by step series of examples that tell you how to get a development env running

```
Before running any programs, make sure you have your paths (to the code and executables) set up correctly. This is done from one 
location: MutateCompute/configuration/pathsAndDictionaries.yaml. Here just add new absolute paths or edit the following 
existing ones: 
Four executables on cluster:

> path_zeus_r_exe
> path_zeus_fx_exe
> path_zeus_agad_exe
> path_zeus_qsub_exe

Three executables on local machine: 
> path_local_r_exe
> path_local_fx_exe
> path_local_agad_exe

One path to project on cluster:
> path_zeus_snpeffect
One path to project on local machine:
> path_local_pyproj_mutatecompute
And any others you might want, such as a path to pdb or fasta files if these are kept outside of the project directory such as:
> path_local_pdb_fasta_repository

```

```
To run FoldX BuildModel mutation, go to src/launchers/KickOffManual.py file and select the pdb file(s) to run, set 
'do_foldx_buildmodel' value to True in 'operations' dictionary. Select the amino acids you want to mutate your protein to. 
For example, to use all 20 amino acids reference a list from an enum as follows: enter the following import line at top of script 
from src.enums.AminoAcids import AA. Now you can access any values from the enum as follows: AA.LIST_ALL_20_AA.value.
  
```


```
If running locally, you can start the python launcher script directly from the IDE or via sh KickOffManual.py in terminal. The 
paths should automatically be set to local paths. 
If running remotely (on cluster), start using sh KickOffZeus.sh which is in switchlab/group/shazib/SnpEffect/bash directory.

(Make sure you have copied over any changes to the cluster using for example:
scp -r -P 7788 /Users/u0120577/PycharmProjects/MutateCompute/src shazib@zeus.psb.ugent.be:/switchlab/group/shazib/SnpEffect to 
copy over all contents of the src/ directory over. This overwrites other modules and files with the same name in the same 
location but will not delete any other older files with different names of in different locations.) 

```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Deployment is not applicable. 

## Built With

No build management tools were used. No web frameworks were used. 

## Contributing

I currently am not expecting any outside contributors for this project in its current form. However if this changes in the 
future I would likely refer to others pre-written guidelines such as read [CONTRIBUTING.md] (https://gist.github
.com/PurpleBooth/b24679402957c63ec426) for details on code of conduct, and the process for submitting pull requests.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Shahin Zibaee** - *Informed and adapted from Python scripts written by Rob van der Kant* - [originalScriptsAtJan2018]
(https://github.com/snphive/originalScriptsAtJan2018)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Rob van der Kant
* Mark Viers (help with setting up remote debugging) 