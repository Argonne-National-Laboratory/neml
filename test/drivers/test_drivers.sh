#!/bin/sh

ifile="../regression/models.txt"

IWD=$(pwd)
cd "$(dirname $0)"

while IFS= read -r line
do
      echo "Beginning test: $line"
      printf "\tTesting cxx driver..."
      ../../util/cxx_interface/cxxsimple ../regression/reference.xml "$line" 0.02 100.0 20 300
      printf " done\n"
      printf "\tTesting c driver..."
      ../../util/c_interface/csimple ../regression/reference.xml "$line" 0.02 100.0 20 300
      printf " done\n"
      printf "\tTesting c driver..."
      ../../util/f_interface/fsimple ../regression/reference.xml "$line" 0.02 100.0 20 300
      printf " done\n"
      echo ""
done < "$ifile"

cd $IWD
