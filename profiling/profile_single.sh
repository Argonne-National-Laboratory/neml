#!/bin/bash

ifile="../test/regression/models.txt"

if [[ $# -ne 1 ]]; then
    echo "Script takes a single parameter, the model name"
    echo "Options are:"
    while IFS= read -r line
    do
      echo "    $line"
    done < "$ifile"
    exit 2
fi

valgrind --tool=callgrind ../util/cxx_interface/cxxsimple ../test/regression/reference.xml "$1" 0.02 100.0 20 300
