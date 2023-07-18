

FAZ: A flexible auto-tuned modular error-bounded compression framework for scientific data

## Introduction

This is the source code of FAZ: A flexible auto-tuned modular error-bounded compression framework for scientific data

## Dependencies

Please Installing the following dependencies:

* Python >= 3.6
* numpy 
* PyWavelets

## 3rd party libraries/tools

* Zstandard (https://facebook.github.io/zstd/). Not mandatory to be mannually installed as Zstandard v1.4.5 is included and will be used if libzstd can not be found by
  pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include

##Command line executable Usage

You can use the executable 'faz' command to do the compression/decompression. Just run "faz" command without any argument to check the instructions for its usage.
faz command is quite like sz3 command, please add the -F argument to activate the FAZ features with default hyper-parameters or customize those hyper-parameters via a config file (not adding -F will make faz fall back to sz3).






