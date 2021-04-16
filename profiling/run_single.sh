#!/bin/bash

ifile="models.txt"

if [[ $# -ne 1 ]]; then
    echo "Script takes a single parameter, the model name"
    echo "Options are:"
    while IFS= read -r line
    do
      echo "    $line"
    done < "$ifile"
    exit 2
fi

../util/cxx_interface/cxxsimple reference.xml "$1" 0.04 100.0 100 300
