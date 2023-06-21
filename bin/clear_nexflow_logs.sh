#!/bin/bash

# Delete all .nextflow.log files
rm -f .nextflow.log*

# Call ./clear_nexflow_logs.sh "-all" to clean data folders as well
if [[ $# -gt 0 && "$1" == "-all" ]]; then
    rm -rf .nextflow/*
    rm -rf work/*
fi
