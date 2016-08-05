#! /usr/bin/env bash

# requires googler-perftools
# "gperftools"
export LD_PRELOAD=/usr/local/lib/libprofiler.dylib
CPUPROFILE="myprof.log" R -f standard-run.R
