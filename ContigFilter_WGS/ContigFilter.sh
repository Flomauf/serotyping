#!/bin/bash

if [ ! -d "Output" ]; then
mkdir Output
fi

echo 'Running FileConverter.py'
python bin/FileConverter.py

echo 'Running AssemblyStat.py'
python bin/AssemblyStat.py

echo 'Running Quast.py'
python bin/Quast.py > /dev/null
