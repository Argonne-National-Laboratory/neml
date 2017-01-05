#!/bin/bash
FC=ifort
INCDIR=../../src
LDIR=../../lib

$FC -O0 -g -I $INCDIR -L $LDIR -lneml -Wl,-rpath=$LDIR fsimple.f -o fsimple
