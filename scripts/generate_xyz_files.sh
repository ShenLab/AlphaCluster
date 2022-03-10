#!/bin/bash

for filename in ../proteins/*.pdb.gz
do
    echo $(basename $filename .pdb.gz)
    gunzip -k $filename
    python xplo2xyz.py "../proteins/$(basename $filename .pdb.gz).pdb" "../proteins/$(basename $filename .pdb.gz).xyz"
    rm "../proteins/$(basename $filename .pdb.gz).pdb"
done
