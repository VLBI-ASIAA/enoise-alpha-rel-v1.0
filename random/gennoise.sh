#!/bin/bash
if [ ! -d "ndata" ]; then
  mkdir ndata
fi

cd ndata

$ENOISE_HOME/random/genran \
  -s 58325 -h 5 \
  -n common.dat:2.0 \
  -n delay.dat:0.5 \
  -n noise0.dat:2.0 \
  -n noise1.dat:2.0 \
  -n noise2.dat:2.0
