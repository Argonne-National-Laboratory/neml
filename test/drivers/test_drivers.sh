#!/bin/sh

cd "$(dirname $0)"

for d in ../test_regression/test_*/
do
	echo "Beginning $(basename $d)"
	printf "\tTesting cxx driver..."
	if ! ../../util/cxx_interface/cxxsimple "$d/model.xml" "model" 0.1 100.0 100 300; then
		printf " error!\n"
		exit -1
	fi
	printf " done\n"
	printf "\tTesting c driver..."
	if ! ../../util/c_interface/csimple "$d/model.xml" "model" 0.1 100.0 100 300; then
		printf " error!\n"
		exit -1
	fi
	printf " done\n"
	printf "\tTesting fortran driver..."
	if ! ../../util/f_interface/fsimple "$d/model.xml" "model" 0.1 100.0 100 300; then
		printf " error!\n"
		exit -1
	fi
	printf " done\n"
	echo ""
done
