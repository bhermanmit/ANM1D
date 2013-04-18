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

    def extract_mean(self,tally_id,score_id,nuclide_id=-1):
        # extract results
        results = self.sp.extract_results(tally_id,score_id,nuclide_id)

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
        if surf == 'left':
            self.currL = np.zeros(self.ng)
        else:
            self.currR = np.zeros(self.ng)
        partial_currents = self.extract_mean(tally_id,score_id)
        for i in range(self.ng):
            if surf == 'left':
                self.currL[i] = partial_currents[1,i,1,1,0] - \
                               partial_currents[0,i,1,1,0]
                print partial_currents[1,i,1,1,0],partial_currents[0,i,1,1,0]
                print 'Left Current:',self.currL[i]
            elif surf == 'right':
                self.currR[i] = partial_currents[1,i,1,1,1] - \
                               partial_currents[0,i,1,1,1]
                print partial_currents[1,i,1,1,1],partial_currents[0,i,1,1,1]
                print 'Right Current:',self.currR[i]
        if surf == 'left':
            self.currL = self.currL[::-1]
        else:
            self.currR = self.currR[::-1]

    def write_matlab_binary(self):
        mdict = { 'flux':self.flux,
                  'meshflux':self.meshflux,
                  'meshnsf':self.meshnufission,
                  'sigt':self.totxs,
                  'sigs':self.scat,
                  'nsigf':self.nfissvec,
                  'chi':self.chi,
                  'diff':self.diff,
                  'currL':self.currL,
                  'currR':self.currR,
                  'keff':self.keff,
                  'H1totalrate':self.H1totalrate,
                  'totalrate':self.totalrate,
                  'H1p1rate':self.H1p1rate,
                  'p1rate':self.p1rate,
                  'fluxrate':self.fluxrate }
        scipy.io.savemat(self.name+'.mat',mdict,appendmat=False)

    def process_meshflux(self,tally_id,score_id):
        self.meshflux = self.extract_mean(tally_id,score_id)
        self.meshflux = self.meshflux[:,:,0,0]
        for i in range(np.size(self.meshflux,1)):
            self.meshflux[:,i] = self.meshflux[::-1,i]

    def process_meshnufission(self,tally_id,score_id):
        self.meshnufission = self.extract_mean(tally_id,score_id)
        self.meshnufission = self.meshnufission[:,:,0,0]
        for i in range(np.size(self.meshnufission,1)):
            self.meshnufission[:,i] = self.meshnufission[::-1,i]

    def process_H1totalrate(self,tally_id,score_id,nuclide_id=-1):
        self.H1totalrate = self.extract_mean(tally_id,score_id,nuclide_id)
        self.H1totalrate = self.H1totalrate[:,0,0,0]
        self.H1totalrate = self.H1totalrate[::-1]

    def process_totalrate(self,tally_id,score_id,nuclide_id=-1):
        self.totalrate = self.extract_mean(tally_id,score_id,nuclide_id)
        self.totalrate = self.totalrate[:,0,0,0]
        self.totalrate = self.totalrate[::-1]

    def process_H1p1rate(self,tally_id,score_id,nuclide_id=-1):
        self.H1p1rate = self.extract_mean(tally_id,score_id,nuclide_id)
        self.H1p1rate = self.H1p1rate[:,0,0,0]
        self.H1p1rate = self.H1p1rate[::-1]

    def process_p1rate(self,tally_id,score_id,nuclide_id=-1):
        self.p1rate = self.extract_mean(tally_id,score_id,nuclide_id)
        self.p1rate = self.p1rate[:,0,0,0]
        self.p1rate = self.p1rate[::-1]

    def process_fluxrate(self,tally_id,score_id,nuclide_id=-1):
        self.fluxrate = self.extract_mean(tally_id,score_id,nuclide_id)
        self.fluxrate = self.fluxrate[:,0,0,0]
        self.fluxrate = self.fluxrate[::-1]


def main(statepoint_file_both, statepoint_file_SA):

    # read in statepoint file
    sp_both = statepoint.StatePoint(statepoint_file_both)
    sp_both.read_results()
    sp_SA = statepoint.StatePoint(statepoint_file_SA)
    sp_SA.read_results()

    # create mesh data objects
    reg1 = MeshData(sp_SA, 'reg1')
    reg2 = MeshData(sp_both, 'reg2')

    # read in left mesh tallies
    reg1.process_flux(1,'flux')
    reg1.process_totalxs(1,'total')
    reg1.process_diff(1,'scatter-n')
    reg1.process_nuscattmat(2,'nu-scatter')
    reg1.process_nufissmat(2,'nu-fission')
    reg1.calc_nfissvec()
    reg1.currL = [0.0,0.0]
    reg1.currR = [0.0,0.0]
    reg1.process_meshflux(5,'flux')
    reg1.process_meshnufission(5,'nu-fission')
    reg1.process_H1totalrate(3,'total',1001)
    reg1.process_totalrate(3,'total')
    reg1.process_H1p1rate(3,'scatter-n',1001)
    reg1.process_p1rate(3,'scatter-n')
    reg1.process_fluxrate(4,'flux')

    # read in right mesh tallies
    reg2.process_flux(2,'flux')
    reg2.process_totalxs(2,'total')
    reg2.process_diff(2,'scatter-n')
    reg2.process_nuscattmat(4,'nu-scatter')
    reg2.process_nufissmat(4,'nu-fission')
    reg2.calc_nfissvec()
    reg2.process_current(6,'current','left')
    reg2.process_current(6,'current','right')
    reg2.process_meshflux(11,'flux')
    reg2.process_meshnufission(11,'nu-fission')
    reg2.process_H1totalrate(8,'total',1001)
    reg2.process_totalrate(8,'total')
    reg2.process_H1p1rate(8,'scatter-n',1001)
    reg2.process_p1rate(8,'scatter-n')
    reg2.process_fluxrate(10,'flux')

    # calculate keff from tallies
    reg1.keff = 1.21293
    reg2.keff = 1.02926

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
    statepoint_file_both = sys.argv[1]
    statepoint_file_SA = sys.argv[2]
    main(statepoint_file_both, statepoint_file_SA)
