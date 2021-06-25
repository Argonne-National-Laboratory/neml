#!/bin/sh

cd "$(dirname $0)"

for d in ../test_regression/test_*/
do
	echo "Beginning $(basename $d)"
	printf "\tTesting cxx driver..."
	../../util/cxx_interface/cxxsimple "$d/model.xml" "model" 0.1 100.0 100 300
	printf " done\n"
	printf "\tTesting c driver..."
	../../util/c_interface/csimple "$d/model.xml" "model" 0.1 100.0 100 300
	printf " done\n"
	printf "\tTesting fortran driver..."
	../../util/f_interface/fsimple "$d/model.xml" "model" 0.1 100.0 100 300
	printf " done\n"
	echo ""
done
