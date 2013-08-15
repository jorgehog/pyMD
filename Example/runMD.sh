#!/usr/bin/env bash

cp pyMD_main.py pyMD_main.pyx
python setup.py build_ext --inplace
rm pyMD_main.c 
rm pyMD_main.pyx
python -c "import pyMD_main"

