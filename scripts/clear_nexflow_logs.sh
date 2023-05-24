#!/bin/bash

# Delete all .nextflow.log files
rm -f .nextflow.log*

# Call ./clear_nexflow_logs.sh "mode=all" to clean data folders as well
if [[ $# -gt 0 && "$1" == "mode=all" ]]; then
    rm -rf .nextflow/*
    rm -rf work/*
fi
