#!/bin/bash

DLL_PATH=$1

if [[ -z ${DLL_PATH} ]]; then
    echo -e "Missing CTDSE DLL path as argument.\nCall ./run_tests.sh <path-to-CTDSE-DLL>"
    echo "Stopping."
    exit 1
fi

python3 tests/tests.py -d ${DLL_PATH}