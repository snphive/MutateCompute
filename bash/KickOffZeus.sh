#!/usr/bin/env bash
source /home/shazib/.bash_profile
source /switchlab/group/shazib/SnpEffect/venv/bin/activate

export PYTHONPATH=$PYTHONPATH:$PWD/src/pycharm-debug.egg

#!/usr/bin/env python3
python3 /switchlab/group/shazib/SnpEffect/src/KickOff.py 'use_cluster=True'