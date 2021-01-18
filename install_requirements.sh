#!/bin/bash

conda config --add channels conda-forge

while read requirement; do
    echo "$requirement"
    conda install --yes $requirement
done < requirements.txt
