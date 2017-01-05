#!/bin/bash
CC=gcc
INCDIR=../../src
LDIR=../../lib

$CC -O0 -g -I $INCDIR -L $LDIR -lneml -Wl,-rpath=$LDIR csimple.c -o csimple
