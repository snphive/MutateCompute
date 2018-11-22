#!/usr/bin/env bash
source /home/shazib/.bash_profile
source /switchlab/group/shazib/SnpEffect/venv/bin/activate

export PYTHONPATH=$PYTHONPATH:/switchlab/group/shazib/SnpEffect/venv/pycharm-debug.egg
export PYTHONPATH=$PYTHONPATH:/switchlab/group/shazib/SnpEffect/venv/lib/python3.7/site-packages
export PYTHONPATH=$PYTHONPATH:/switchlab/group/shazib/SnpEffect

#!/usr/bin/env python3
python3 /switchlab/group/shazib/SnpEffect/src/launchers/BlastOff.py 'use_cluster=True'

