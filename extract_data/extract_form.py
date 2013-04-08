#!/usr/bin/python
# Filename: extract_data.py

# import packages
import statepoint
import numpy as np
import scipy.io
import os
import sys

class MeshData:
    def __init__(self,sp,name):
        self.sp = sp
        self.name = name

    def extract_mean(self,tally_id,score_id):
        # extract results
        results = self.sp.extract_results(tally_id,score_id)

        # copy means over
        mean = results['mean'].copy()

        # reshape the mean
        mean = mean.reshape(results['bin_max'],order='F')

        return mean

    def process_flux(self,tally_id,score_id):
        self.flux = self.extract_mean(tally_id,score_id)
        self.flux = self.flux[:,:,0,0]
        for i in range(np.size(self.flux,1)):
          self.flux[:,i] = self.flux[::-1,i]

    def write_matlab_binary(self):
        mdict = { 'flux':self.flux }
        scipy.io.savemat(self.name+'.mat',mdict,appendmat=False)


def main(statepoint_file_low, statepoint_file_high):

    # read in statepoint file
    sp_low = statepoint.StatePoint(statepoint_file_low)
    sp_high = statepoint.StatePoint(statepoint_file_high)
    sp_low.read_results()
    sp_high.read_results()

    # create mesh data objects
    reg1 = MeshData(sp_low, 'reg1form')
    reg2 = MeshData(sp_high, 'reg2form')

    # read in left mesh tallies
    reg1.process_flux(1,'flux')

    # read in right mesh tallies
    reg2.process_flux(1,'flux')

    # write matlab output binaries
    reg1.write_matlab_binary()
    reg2.write_matlab_binary()

if __name__ == "__main__":
    statepoint_file_low = sys.argv[1]
    statepoint_file_high = sys.argv[2]
    main(statepoint_file_low, statepoint_file_high)
