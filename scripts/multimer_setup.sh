#!/bin/bash

ID=$1
wget https://files.rcsb.org/download/$ID.pdb -P multimers/
python scripts/xplo2xyz_multimer.py multimers/$ID.pdb 
