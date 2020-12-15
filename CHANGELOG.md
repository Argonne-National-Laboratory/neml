# Change Log

## 1.3.0 - 9/15/2020
* Removed boost dependency
* Added kinematic hardening to crystal plasticity model
* Added work-based damage model
* Block conversion routines for going to/from tensor notation
* Updated build/release process

## 1.2.2 - 4/7/2020
* Bugfix for new python package distribution

## 1.2.1 - 4/7/2020
### Bugfix
* Documentation issue between breathe and sphinx versions resolved

## 1.2.0 - 4/7/2020
### Features
* Added an entire crystal plasticity subsystem (submodule cp) for single crystal constitutive response
* Added a python package distribution with a limited substep of the library featuers
* Test suite now includes regression tests, in addition to the unit tests
* Basic Windows compilation instructions in the docs and tested
* Simplified high level material model integration algorithm

## 1.1.0 - 4/24/2019
### Features
* Switch to python3
* Added large deformations via a Truesdell objective rate
* Moved from libxml++ to rapidxml to reduce dependencies
