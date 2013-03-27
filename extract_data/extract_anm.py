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
        self.flux = self.flux[:,0,0,0]
        self.ng = np.size(self.flux,0)
        self.flux = self.flux[::-1]

    def process_totalxs(self,tally_id,score_id):
        self.totxs = self.extract_mean(tally_id,score_id)
        self.totxs = self.totxs[:,0,0,0]
        self.totxs = self.totxs[::-1]
        self.totxs = self.totxs/self.flux

    def process_diff(self,tally_id,score_id):
        self.diff = self.extract_mean(tally_id,score_id)
        self.diff = self.diff[:,0,0,0]
        self.diff = self.diff[::-1]
        self.diff = self.diff/self.flux
        self.diff = 1.0/(3*(self.totxs-self.diff))

    def process_nuscattmat(self,tally_id,score_id):
        scat = self.extract_mean(tally_id,score_id)
        ng = self.ng
        self.scat = np.zeros((ng,ng)) 
        for i in range(ng):
            for j in range(ng):
                self.scat[ng-j-1,ng-i-1] = scat[i,j]/self.flux[ng-j-1]

    def process_nufissmat(self,tally_id,score_id):
        nfiss = self.extract_mean(tally_id,score_id)
        ng = self.ng
        self.nfiss = np.zeros((ng,ng))
        for i in range(ng):
            for j in range(ng):
                self.nfiss[ng-j-1,ng-i-1] = nfiss[i,j]/self.flux[ng-j-1]

    def calc_nfissvec(self):
        totnfiss = self.nfiss
        self.nfissvec = np.zeros(self.ng)
        self.chi      = np.zeros(self.ng)
        for i in range(self.ng):
            totnfiss = np.sum(totnfiss,0)
        for i in range(self.ng):
            self.nfissvec[i] = sum(self.nfiss[i,:])
            self.chi[i] = sum(self.nfiss[:,i])/totnfiss

    def process_current(self,tally_id,score_id,surf):
        self.curr = np.zeros(self.ng)
        partial_currents = self.extract_mean(tally_id,score_id)
        for i in range(self.ng):
            if surf == 'left':
                self.curr[i] = partial_currents[1,i,1,1,0] - \
                               partial_currents[0,i,1,1,0]
            elif surf == 'right':
                self.curr[i] = partial_currents[1,i,1,1,1] - \
                               partial_currents[0,i,1,1,1]
        self.curr = self.curr[::-1]

    def write_matlab_binary(self):
        mdict = { 'flux':self.flux,
                  'sigt':self.totxs,
                  'sigs':self.scat,
                  'nsigf':self.nfissvec,
                  'chi':self.chi,
                  'diff':self.diff,
                  'curr':self.curr,
                  'keff':self.keff  }
        scipy.io.savemat(self.name+'.mat',mdict,appendmat=False)


def main(statepoint_id):

    # read in statepoint file
    sp = statepoint.StatePoint('statepoint.'+statepoint_id+'.binary')
    sp.read_results()

    # create mesh data objects
    reg1 = MeshData(sp, 'reg1')
    reg2 = MeshData(sp, 'reg2')

    # read in left mesh tallies
    reg1.process_flux(1,'flux')
    reg1.process_totalxs(1,'total')
    reg1.process_diff(1,'scatter-n')
    reg1.process_nuscattmat(3,'nu-scatter')
    reg1.process_nufissmat(3,'nu-fission')
    reg1.calc_nfissvec()
    reg1.process_current(5,'current','right')

    # read in right mesh tallies
    reg2.process_flux(2,'flux')
    reg2.process_totalxs(2,'total')
    reg2.process_diff(2,'scatter-n')
    reg2.process_nuscattmat(4,'nu-scatter')
    reg2.process_nufissmat(4,'nu-fission')
    reg2.calc_nfissvec()
    reg2.process_current(6,'current','left')

    # calculate keff from tallies
    keff = calc_keff(reg1,reg2)
    reg1.keff = keff
    reg2.keff = keff

    # write matlab output binaries
    reg1.write_matlab_binary()
    reg2.write_matlab_binary()

def calc_keff(reg1,reg2):
    reg1nsigf = np.sum(reg1.nfissvec*reg1.flux)
    reg2nsigf = np.sum(reg2.nfissvec*reg2.flux)
    reg1nabs = np.sum((reg1.totxs - np.sum(reg1.scat,1))*reg1.flux)
    reg2nabs = np.sum((reg2.totxs - np.sum(reg2.scat,1))*reg2.flux)

    keff = (reg1nsigf + reg2nsigf)/(reg1nabs + reg2nabs)

    return keff

if __name__ == "__main__":
    statepoint_id = sys.argv[1]
    main(statepoint_id)
