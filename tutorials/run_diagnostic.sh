#!/bin/bash

# check for python
if command -v python3 >/dev/null 2>&1; then
    PYTHONCMD="python3"
    python3 --version 
elif command -v python >/dev/null 2>&1; then
    PYTHONCMD="python"
    python --version
else
    echo "ERROR: python is not installed"
    exit 1
fi
echo "using $PYTHONCMD"

export PYTHONCMD

# check for python libs 
$PYTHONCMD .pylib_diagnostic.py

# just check for awk... never know
if command -v awk >/dev/null 2>&1; then
    echo "awk is installed"
else
    echo "ERROR: awk is not installed"
    exit 1
fi

