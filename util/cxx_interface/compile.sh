#!/bin/bash
CXX=g++
INC="-I../../src -I/usr/include/libxml++-2.6 -I/usr/lib64/libxml++-2.6/include -I/usr/include/glibmm-2.4 -I/usr/lib64/glibmm-2.4/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -I/usr/include/sigc++-2.0 -I/usr/lib64/sigc++-2.0/include"
LDIR=../../lib

$CXX -O0 -g -std=gnu++11 $INC -L $LDIR -lneml -Wl,-rpath=$LDIR cxxsimple.cxx -o cxxsimple
