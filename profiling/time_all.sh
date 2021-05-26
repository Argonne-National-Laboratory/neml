#!/bin/sh

ifile="models.txt"
repeats=10

i=0
while [ $i -lt $repeats ]
do
      while IFS= read -r line
      do
            ../util/cxx_interface/cxxsimple reference.xml "$line" 0.04 100.0 100 300
      done < "$ifile"
      i=$(( $i + 1 ))
done
