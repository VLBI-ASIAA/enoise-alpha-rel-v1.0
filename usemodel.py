#!/usr/bin/python

# This script is a modification of the original code
# provided by Adam Deller

import sys, os
import difxmodel

model = difxmodel.Model(sys.argv[1])
mjd = model.startmjd
seconds = model.startsec
sourcename = model.sources[0].name

#for num in range(len(model.antennas)):
#  for sec in range(int(sys.argv[2])):
#    csec = seconds + sec
#    delay, rate = model.getDelayAndRate(mjd, csec, model.antennas[num].name, sourcename)


if not os.path.exists("drate"):
  os.makedirs("drate")
for num in range(len(model.antennas)):
  file = "drate/rate" + str(num) + ".dat"
  f = open(file, 'w')
  
  # read seconds from command line
  for sec in range(int(sys.argv[2])):
    csec = seconds + sec
    delay, rate, accel = model.getDelayAndRate(mjd, csec, model.antennas[num].name, sourcename)
    f.write("%.9f %.9f %.9f\n" % (delay, rate, accel))
  
  f.close()
