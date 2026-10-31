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
      description='Script to print out decays of particle containing specific daughters')
parser.add_argument('-decfile', default='DECAY.DEC')
parser.add_argument('-particle', default='')
parser.add_argument('-daughter1', default=[], nargs='*')
parser.add_argument('-daughter2', default=[], nargs='*')

args = parser.parse_args()

print(args)

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
    # Now check whether we have given daughters
    wanted = False
    if len(args.daughter1)!=0:
      for ii in args.daughter1:
        if ii in line: # one of the daughers is in
          if len(args.daughter2)!=0:
            for jj in args.daughter2:
              if jj in line:
                wanted = True
          else:
            wanted = True
    if wanted:
      print(line.strip())

decFiles.close()
   
