#!/usr/bin/env python

########################################################################
# Copyright 1998-2022 CERN for the benefit of the EvtGen authors       #
#                                                                      #
# This file is part of EvtGen.                                         #
#                                                                      #
# EvtGen is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# EvtGen is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     #
########################################################################

import os,sys
import argparse

parser = argparse.ArgumentParser(
      description='Script to count total BF for given particle')
parser.add_argument('-decfile', default='DECAY.DEC')
parser.add_argument('-particle', default='')

args = parser.parse_args()

decFiles = open( args.decfile )

inDecay = False
bfSum = 0.0
nDecays = 0
for line in decFiles:
  if (inDecay ==  False) and (line.find('Decay') != -1):
    ss = line.split(' ')
    if ( ss[len(ss)-1].strip() == args.particle ):
      inDecay = True
      print(line)
      continue
  if (inDecay == True) and (line.find('Enddecay') != -1):
    inDecay = False
    break
  if inDecay == True:
    if line[0] == '#' or len(line)<2: # This is comment, skip
      continue
    print(line.strip())
    bfSum += float(line.strip().split(' ')[0])
    nDecays += 1

print('Counted ',nDecays,' decays with total BF = ',bfSum)
print('Missing to 1 is ',1-bfSum)
decFiles.close()
   
