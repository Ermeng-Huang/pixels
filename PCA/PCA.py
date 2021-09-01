# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:49:17 2020

@author: hem
"""
import sys
sys.path.insert(1,'d:/code')
import os
import csv
import h5py
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pixelStats import baselineVector
from Arrow3D import Arrow3D
from matplotlib import rcParams
from scipy.spatial.distance import euclidean
from sklearn.model_selection import KFold
import scipy.stats as stats
import re


class PCA_stats:
    def __init__(self):
        self.su_list = []
        self.avail = []
        self.kf = KFold(2, shuffle=True)

    def splitMean(self, FR, su_idx, mm, std, random_type=None):
        if random_type is None:
            out = np.squeeze(np.mean((FR[su_idx,:, :] - mm) / std, axis=1))
        elif random_type == 'half':
            sel = list(self.kf.split(np.arange(FR.shape[1])))[0]
            out = np.squeeze(np.mean((FR[su_idx, :,sel[0],] - mm) / std, axis=1))
        else:
            sel = np.random.choice(FR.shape[2],1)
            out = np.squeeze((FR[su_idx, :, sel, ] - mm) / std)
        return out

        ### unprocessed data entry point
    
                
    def processCTDStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None, random_type=None):
        
        ### TODO: when variables are empty
        if welltrain_window:            
            FR_S1_laseron = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 4, trials[8, :] == 2, welltrain_window, correct_resp,)),axis=0, ), ]
            FR_S2_laseron = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 8, trials[8, :] == 2, welltrain_window, correct_resp,)),axis=0, ), ]
            FR_S1_laseroff = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 4, trials[8, :] == -1, welltrain_window, correct_resp,)),axis=0, ), ]
            FR_S2_laseroff = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 8, trials[8, :] == -1, welltrain_window, correct_resp,)),axis=0, ), ]
        else:
            FR_S1_laseron = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 4, trials[8, :] == 2)),axis=0, ), ]
            FR_S2_laseron = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 8, trials[8, :] == 2)),axis=0, ), ]
            FR_S1_laseroff = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 4, trials[8, :] == -1)),axis=0, ), ]
            FR_S2_laseroff = trial_FR[:, :, np.all(np.vstack((trials[4, :] == 8, trials[8, :] == -1)),axis=0, ), ]
        for su_idx in range(trial_FR.shape[0]):
            (mm, std) = baselineVector(trial_FR[su_idx,:, :])
            if std > 0 and np.all([x.shape[2] >= 2 for x in (
                    FR_S1_laseron,FR_S2_laseron,FR_S1_laseroff,FR_S2_laseroff
            )]):
                self.su_list.append(
                    {
                        "S1_laseron": self.splitMean(FR_S1_laseron, su_idx, mm, std, random_type),
                        "S2_laseron": self.splitMean(FR_S2_laseron, su_idx, mm, std, random_type),
                        "S1_laseroff": self.splitMean(FR_S1_laseroff, su_idx, mm, std, random_type),
                        "S2_laseroff": self.splitMean(FR_S2_laseroff, su_idx, mm, std, random_type),                        
                    }
                )
                self.avail.append(True)
            else:
                self.avail.append(False)

    def get_features(self):
        # breakpoint()
        return (self.su_list, self.avail)


def traverse(path):
    for (basepath, dirs, files) in os.walk(path):
        if "cluster_info.tsv" in files:
            yield basepath
            
def judgePerformance(trials, criteria=75):
    """

    Parameters
    ----------
    trials : TYPE
        behaviral trials array.

    Returns
    -------
    str
        readable description of the learning stage.
    int
        single integer code for the learning stage.
    trials
        if well-trained, the trials in the engaging window, all trials otherwise.

    """
    sample_loc = 4 if trials.shape[0] > 6 else 2
    test_loc = 5 if trials.shape[0] > 6 else 3
    lick_loc = 6 if trials.shape[0] > 6 else 4
    
    if trials.shape[1] >= 40:
        correctResp = np.bitwise_xor(trials[sample_loc, :] == trials[test_loc, :], trials[lick_loc, :] == 1)
        inWindow = np.zeros((trials.shape[1],), dtype="bool")
        i = 40
        while i < trials.shape[1]:
            #            if np.sum(correctResp[i-40:i])>=32:
            if np.sum(correctResp[i - 40: i]) >= criteria * 40 / 100:
                inWindow[i - 40: i] = 1
            i += 1
        if np.sum(inWindow) >= 40:  # Well Trained
            return ("wellTrained", 3, inWindow, correctResp)

        else:
            inWindow = np.zeros((trials.shape[1],), dtype="bool")
            licks = trials[:, lick_loc] == 1
            i = 40
            while i < trials.shape[1]:
                if np.sum(licks[i - 40: i]) >= 16:  # Learning
                    inWindow[i - 40: i] = 1
                i += 1
            if np.sum(inWindow) >= 40:
                return ("learning", 2, inWindow, correctResp)
            elif np.sum(licks) <= trials.shape[1] // 10:  # Passive
                return ("passive", 0, np.ones_like(inWindow), correctResp)
            else:
                return ("transition", 1, np.ones_like(inWindow), correctResp)
            
def get_dataset(denovo=False, to_explore=True, random_type=None):
    features_per_su = []
    avails = []
    reg_list = []
    suid_reg=[];
    SU_id_list=[];
    session=os.listdir('D:\pixel-optogenetic\DataSum')
    if denovo:        
        for path in sorted(traverse('D:\pixel-optogenetic\DataSum')):
            # print(path)  
            rootpath=re.search(r"(?<=DataSum\\)(.*)", path)[1]
            s=float(session.index(re.search(r"(.*)(?=\\)", rootpath)[1]))+1  # session ID
            t=float(re.search(r"(imec)(.*)", rootpath)[2][0]) # track ID in per session
            with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
                # print(list(ffr.keys()))
                if not "SU_id" in ffr.keys():                    
                    print("missing su_id key in path ", path)
                    continue
                dset = ffr["SU_id"]
                SU_ids = np.array(dset, dtype="double")
                dset = ffr["FR_All"]
                trial_FR = np.array(dset, dtype="double")
                dset = ffr["Trials"]
                trials = np.array(dset, dtype="double")            
                
            # (_perf_desc, perf_code, welltrain_window, correct_resp,) = judgePerformance(trials)
            # if perf_code != 3:
            #     continue
            
            with open(os.path.join(path, "su_id2reg_allclass.csv")) as csvfile:
                l = list(csv.reader(csvfile))[1:]  
                suid_reg.extend(l[i] for i in range(len(l)) )
            
            
            suid_reg = []
            with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
                l = list(csv.reader(csvfile))[1:]
                suid_reg = [list(i) for i in zip(*l)]

            currStats = PCA_stats()
            # currStats.processCTDStats(trial_FR, trials, welltrain_window, correct_resp, random_type=random_type)
            currStats.processCTDStats(trial_FR, trials, random_type=random_type)
            (feat, avail) = currStats.get_features()
            features_per_su.extend(feat)
            avails.extend(avail)
            reg_list.extend(suid_reg[1])
            SU_id_list.extend(SU_ids[0]+100000*s+10000*t) 
            ### DEBUG in small subsets
            # if len(features_per_su)>50:
            #     break

        ### save to npz file
        # np.savez_compressed("ctd_stats.npz", features_per_su=features_per_su, reg_list=reg_list, avails=avails)

    ### load from npz file
    else:
        
        fstr = np.load("ctd_stats.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        reg_list = fstr["reg_list"].tolist()
        avails = fstr["avails"].tolist()

    return (features_per_su, SU_id_list, avails)


def plotPCA3D_FR(fr, plot_delay=6, to_plot=True, fig=None, ax=None, the_ref=True, coeff=None):
    
    pcamat = np.hstack(
        [
            np.vstack([x["S1_laseron"][32:40] for x in fr]),
            np.vstack([x["S2_laseron"][32:40] for x in fr]),
            np.vstack([x["S1_laseroff"][32:40] for x in fr]),
            np.vstack([x["S2_laseroff"][32:40] for x in fr]),
        ]
    )
    # np.save(f"pcamat_{delay}_sust.npy", pcamat)
    if the_ref: #
        pca = PCA(n_components=20)
        comp = pca.fit_transform(pcamat.T)
        coeff = pca.components_
        ratio = pca.explained_variance_ratio_
    else:
        pcamat_cent = pcamat - np.expand_dims(np.mean(pcamat, axis=1), axis=1)
        comp = np.matmul(coeff, pcamat_cent).T
    if the_ref:
        alpha = 1
        lw = 1.5
    else:
        alpha = 0.25
        lw = 0.5

    if to_plot:       
    
        for i in range(9, 52): # -1s ~ 10s
            arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                          [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                          lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="b", linestyle=':', alpha=alpha)
            ax.add_artist(arr)
            j = i + 68
            arr = Arrow3D([comp[j - 1, 0], comp[j, 0]], [comp[j - 1, 1], comp[j, 1]],
                          [comp[j - 1, 2], comp[j, 2]], mutation_scale=20,
                          lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="b", linestyle='--', alpha=alpha)
            ax.add_artist(arr)
            m = j + 68
            arr = Arrow3D([comp[m - 1, 0], comp[m, 0]], [comp[m - 1, 1], comp[m, 1]],
                          [comp[m - 1, 2], comp[m, 2]], mutation_scale=20,
                          lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="k", linestyle=':', alpha=alpha)
            ax.add_artist(arr)
            n = m + 68
            arr = Arrow3D([comp[n - 1, 0], comp[n, 0]], [comp[n - 1, 1], comp[n, 1]],
                          [comp[n - 1, 2], comp[n, 2]], mutation_scale=20,
                          lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="k", linestyle='--', alpha=alpha)
            ax.add_artist(arr)
            
        if the_ref:
            for bidx in [0, 68, 136, 204]:
                ax.plot3D(
                    comp[bidx:bidx + 4, 0], comp[bidx:bidx + 4, 1],
                    comp[bidx:bidx + 4, 2], "k.", markersize=1, alpha=0.25
                )  # S1S1

                ax.plot3D(
                    comp[bidx + 4:bidx + 8, 0], comp[bidx + 4:bidx + 8, 1],
                    comp[bidx + 4:bidx + 8, 2], "k+", markersize=4, alpha=0.25
                )
                
    ax.set_xlim(min(comp[:,0]), max(comp[:,0]))
    ax.set_ylim(min(comp[:,0]), max(comp[:,0]))
    ax.set_zlim(min(comp[:,0]), max(comp[:,0]))
    
    
        
        
    return (comp, coeff)


def plot3d_all(plot_delay=3, repeats=5, to_plot=False):
    print('Calculated PCA')
    
    fig_all = plt.figure()
    ax_all = fig_all.add_subplot(111, projection='3d')
    fig_sust = plt.figure()
    ax_sust = fig_sust.add_subplot(111, projection='3d')
    fig_trans = plt.figure()
    ax_trans = fig_trans.add_subplot(111, projection='3d')
    fig_sel = plt.figure()
    ax_sel = fig_sel.add_subplot(111, projection='3d')
    fig_nonsel = plt.figure()
    ax_nonsel = fig_nonsel.add_subplot(111, projection='3d')    
    
   
    with h5py.File("D:\pixel-optogenetic\Selectivity_AIopto_0419.hdf5", "r") as ffr:
            dset=ffr["sus_trans_noPermutaion"]
            sus_trans=np.array(dset, dtype="double")                  
    with h5py.File("D:\pixel-optogenetic\LaserModulation.hdf5", "r") as ffr:
            dset=ffr["rank"]
            LaserModulation=np.array(dset, dtype="double").T
    with h5py.File("D:\pixel-optogenetic\Performance_allneuron_0308.hdf5", "r") as ffr:
            dset=ffr["perf"]
            Perf=np.array(dset, dtype="double").T
    
    # with h5py.File("D:\pixel-optogenetic\Selectivity_AIopto_WT.hdf5", "r") as ffr:
    #     dset=ffr["sus_trans"]
    #     sus_trans=np.array(dset, dtype="double").T        
    #     dset=ffr["rank"]
    #     LaserModulation=np.array(dset, dtype="double").T           
        
    Laser_avail=((LaserModulation[:,16]==1)|(LaserModulation[:,17]==1)|(LaserModulation[:,18]==1)|(LaserModulation[:,19]==1))&(Perf[:,1]==1)
    # Laser_avail=((LaserModulation[:,16]==1))&(Perf[:,1]==1)
    (features_per_su, SU_id_list, avails) = get_dataset(denovo=True, random_type=None)    
    trans_avail = (sus_trans[:,1]==1 & np.array([avails])).T
    sust_avail = (sus_trans[:,0]==1 & np.array([avails])).T
    sel_avail = (((sus_trans[:,0]==1) | (sus_trans[:,1]==1)) & Laser_avail & np.array([avails])).T
    nonsel_avail = ((sus_trans[:,0]!=1) & (sus_trans[:,1]!=1) & (sus_trans[:,2]!=1) & (sus_trans[:,3]!=1) & Laser_avail & np.array([avails])).T
    
    fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]   
    fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]
    fr_sel = [f for (f, t) in zip(features_per_su, sel_avail) if t]   
    fr_nonsel = [f for (f, t) in zip(features_per_su, nonsel_avail) if t]
    # reg_arr = trans_fstr['reg_arr']

    (comp_all, coeff_all) = plotPCA3D_FR(features_per_su, plot_delay, to_plot=False, fig=fig_all, ax=ax_all, the_ref=True)
    (comp_sust, coeff_sust) = plotPCA3D_FR(fr_sust, plot_delay, to_plot=False, fig=fig_sust, ax=ax_sust, the_ref=True)
    (comp_trans, coeff_trans) = plotPCA3D_FR(fr_trans, plot_delay, to_plot=False, fig=fig_trans, ax=ax_trans, the_ref=True)
    (comp_sel, coeff_sel) = plotPCA3D_FR(fr_sel, plot_delay, to_plot=False, fig=fig_sel, ax=ax_sel, the_ref=True)
    (comp_nonsel, coeff_nonsel) = plotPCA3D_FR(fr_nonsel, plot_delay, to_plot=False, fig=fig_nonsel, ax=ax_nonsel, the_ref=True)
     
    comp_all_repeat=np.empty((repeats,1),dtype=object)
    comp_sust_repeat=np.empty((repeats,1),dtype=object)
    comp_trans_repeat=np.empty((repeats,1),dtype=object)
    comp_sel_repeat=np.empty((repeats,1),dtype=object)
    comp_nonsel_repeat=np.empty((repeats,1),dtype=object)    
    for r in range(repeats):
        print(r)
        (features_per_su, reg_list, _avails) = get_dataset(denovo=True, random_type='one_trial')
             
        fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]   
        fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]
        fr_sel = [f for (f, t) in zip(features_per_su, sel_avail) if t]   
        fr_nonsel = [f for (f, t) in zip(features_per_su, nonsel_avail) if t]
        
        (comp_all_repeat0, _c) = plotPCA3D_FR(features_per_su, plot_delay, to_plot=False, fig=fig_all, ax=ax_all, the_ref=False, coeff=coeff_all)
        (comp_sust_repeat0, _c) = plotPCA3D_FR(fr_sust, plot_delay, to_plot=False, fig=fig_sust, ax=ax_sust, the_ref=False, coeff=coeff_sust)
        (comp_trans_repeat0, _c) = plotPCA3D_FR(fr_trans, plot_delay, to_plot=False, fig=fig_trans, ax=ax_trans, the_ref=False, coeff=coeff_trans)
        (comp_sel_repeat0, _c) = plotPCA3D_FR(fr_sel, plot_delay, to_plot=False, fig=fig_sel, ax=ax_sel, the_ref=False, coeff=coeff_sel)
        (comp_nonsel_repeat0, _c) = plotPCA3D_FR(fr_nonsel, plot_delay, to_plot=False, fig=fig_nonsel, ax=ax_nonsel, the_ref=False, coeff=coeff_nonsel)
        
        comp_all_repeat[r,0]=comp_all_repeat0
        comp_sust_repeat[r,0]=comp_sust_repeat0
        comp_trans_repeat[r,0]=comp_trans_repeat0
        comp_sel_repeat[r,0]=comp_sel_repeat0
        comp_nonsel_repeat[r,0]=comp_nonsel_repeat0       
          
        
    np.savez(f'D:\pixel-optogenetic\PCA\late2s\PCA_comp_{plot_delay}.npz', comp_all=comp_all, comp_sust=comp_sust, comp_trans=comp_trans, comp_sel=comp_sel, comp_nonsel=comp_nonsel,
             comp_all_repeat=comp_all_repeat,comp_sust_repeat=comp_sust_repeat,comp_trans_repeat=comp_trans_repeat,
             comp_sel_repeat=comp_sel_repeat,comp_nonsel_repeat=comp_nonsel_repeat)
    with h5py.File('D:\pixel-optogenetic\PCA\late2s\PCA_comp.hdf5', "w") as fw:
        fw.create_dataset('comp_all', data=comp_all.astype('double'))   
        fw.create_dataset('comp_sust', data=comp_sust.astype('double'))       
        fw.create_dataset('comp_trans', data=comp_trans.astype('double')) 
        fw.create_dataset('comp_sel', data=comp_sel.astype('double')) 
        fw.create_dataset('comp_nonsel', data=comp_nonsel.astype('double')) 
        fw.create_dataset('comp_all_repeat', data=comp_all_repeat.tolist()) 
        fw.create_dataset('comp_sust_repeat', data=comp_sust_repeat.tolist()) 
        fw.create_dataset('comp_trans_repeat', data=comp_trans_repeat.tolist()) 
        fw.create_dataset('comp_sel_repeat', data=comp_sel_repeat.tolist()) 
        fw.create_dataset('comp_nonsel_repeat', data=comp_nonsel_repeat.tolist())
    
 
    for oneax in [ax_all, ax_sust, ax_trans, ax_sel, ax_nonsel]:
        oneax.set_xlabel('PC1')
        oneax.set_ylabel('PC2')
        oneax.set_zlabel('PC3')

    for onefig in [fig_all, fig_sust, fig_trans, fig_sel, fig_nonsel ]:
        onefig.set_size_inches(65 / 25.4, 65 / 25.4)
    if to_plot:
        fig_all.savefig(f'traj_all_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        fig_sust.savefig(f'traj_sust_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        fig_trans.savefig(f'traj_trans_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        fig_sel.savefig(f'traj_sel_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        fig_nonsel.savefig(f'traj_nonsel_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        plt.show()    


class Dist_stats:
    def __init__(self, ref_comp):
        self.ref_comp = ref_comp
        self.laseron_dist = [] #on S1-S2
        self.laseroff_dist = []
        self.S1_dist=[] # S1 on-off
        self.S2_dist=[]
        self.ref_dist_S1_laseron = []
        self.ref_dist_S2_laseron= []
        self.ref_dist_S1_laseroff = []
        self.ref_dist_S2_laseroff = []

        self.coeff = []

    def append_data(self, fr, coeff):
        pcamat = np.hstack(
            [
                np.vstack([x["S1_laseron"][32:40] for x in fr]),
                np.vstack([x["S2_laseron"][32:40] for x in fr]),
                np.vstack([x["S1_laseroff"][32:40] for x in fr]),
                np.vstack([x["S2_laseroff"][32:40] for x in fr]),
            ]
        )
        # pca = PCA(n_components=20)
        # comp = pca.fit_transform(pcamat.T)
        # coeff = pca.components_
        # ratio = pca.explained_variance_ratio_

        pcamat_cent = pcamat - np.expand_dims(np.mean(pcamat, axis=1), axis=1)
        comp = np.matmul(coeff, pcamat_cent).T

        laseron = []
        laseroff = []
        S1 = []
        S2 = []        
        ref_S1_laseron = []
        ref_S2_laseron = []          
        ref_S1_laseroff = []
        ref_S2_laseroff = []
        for bin in range(0, int(len(comp)/4)):  # time range
            #laser on
            laseron.append(euclidean(comp[bin, :3], comp[bin + int(len(comp)/4), :3]))  # comp range
            laseroff.append(euclidean(comp[bin + 2*int(len(comp)/4), :3], comp[bin + 3*int(len(comp)/4), :3]))  # comp range
            S1.append(euclidean(comp[bin, :3], comp[bin + 2*int(len(comp)/4), :3]))
            S2.append(euclidean(comp[bin + int(len(comp)/4), :3], comp[bin + 3*int(len(comp)/4), :3]))
            ref_S1_laseron.append(euclidean(comp[bin, :3], self.ref_comp[bin, :3]))
            ref_S2_laseron.append(euclidean(comp[bin + int(len(comp)/4), :3], self.ref_comp[bin + int(len(comp)/4), :3]))            
            ref_S1_laseroff.append(euclidean(comp[bin + 2*int(len(comp)/4) , :3], self.ref_comp[bin +2*int(len(comp)/4), :3]))
            ref_S2_laseroff.append(euclidean(comp[bin + 3*int(len(comp)/4), :3], self.ref_comp[bin + 3*int(len(comp)/4), :3]))
        
        self.laseron_dist.append(laseron)
        self.laseroff_dist.append(laseroff)
        self.S1_dist.append(S1)
        self.S2_dist.append(S2)
        self.ref_dist_S1_laseron.append(ref_S1_laseron)
        self.ref_dist_S2_laseron.append(ref_S2_laseron)    
        self.ref_dist_S1_laseroff.append(ref_S1_laseroff)
        self.ref_dist_S2_laseroff.append(ref_S2_laseroff)

    def get_coeff(self, fr):
        pcamat = np.hstack(
            [
                np.vstack([x["S1_laseron"][32:40] for x in fr]),
                np.vstack([x["S2_laseron"][32:40] for x in fr]),
                np.vstack([x["S1_laseroff"][32:40] for x in fr]),
                np.vstack([x["S2_laseroff"][32:40] for x in fr]),
            ]
        )
        pca = PCA(n_components=20)
        comp = pca.fit_transform(pcamat.T)
        coeff = pca.components_
        # ratio = pca.explained_variance_ratio_

        return coeff

    def get_data(self):
        return (self.laseron_dist, self.laseroff_dist, self.S1_dist, self.S2_dist,
                self.ref_dist_S1_laseron, self.ref_dist_S2_laseron, self.ref_dist_S1_laseroff, self.ref_dist_S2_laseroff)


def dist_all(repeats=2):
    print('Calculated distance')
    
    fstr = np.load('D:\pixel-optogenetic\PCA\late2s\PCA_comp_6.npz')
    comp_all = fstr['comp_all']
    comp_sust = fstr['comp_sust']
    comp_trans = fstr['comp_trans']
    comp_sel = fstr['comp_sel']
    comp_nonsel = fstr['comp_nonsel']
    
    dist_all = Dist_stats(comp_all)
    dist_sust = Dist_stats(comp_sust)
    dist_trans = Dist_stats(comp_trans)
    dist_sel = Dist_stats(comp_sel)
    dist_nonsel = Dist_stats(comp_nonsel)

    (features_per_su, reg_list, avails) = get_dataset(denovo=True, to_explore=False, random_type=None)
    with h5py.File("D:\pixel-optogenetic\Selectivity_AIopto_0419.hdf5", "r") as ffr:
            dset=ffr["sus_trans_noPermutaion"]
            sus_trans=np.array(dset, dtype="double")                   
    with h5py.File("D:\pixel-optogenetic\LaserModulation.hdf5", "r") as ffr:
            dset=ffr["rank"]
            LaserModulation=np.array(dset, dtype="double").T
    with h5py.File("D:\pixel-optogenetic\Performance_allneuron_0308.hdf5", "r") as ffr:
            dset=ffr["perf"]
            Perf=np.array(dset, dtype="double").T
            
    # with h5py.File("D:\pixel-optogenetic\Selectivity_AIopto_WT.hdf5", "r") as ffr:
    #     dset=ffr["sus_trans"]
    #     sus_trans=np.array(dset, dtype="double").T        
    #     dset=ffr["rank"]
    #     LaserModulation=np.array(dset, dtype="double").T    
    Laser_avail=((LaserModulation[:,16]==1)|(LaserModulation[:,17]==1)|(LaserModulation[:,18]==1)|(LaserModulation[:,19]==1))& (Perf[:,1]==1)
        
    trans_avail = (sus_trans[:,1]==1 & np.array([avails])).T
    sust_avail = (sus_trans[:,0]==1 & np.array([avails])).T
    sel_avail = (((sus_trans[:,0]==1) | (sus_trans[:,1]==1)) & Laser_avail & np.array([avails])).T
    nonsel_avail = ((sus_trans[:,0]!=1) & (sus_trans[:,1]!=1) & (sus_trans[:,2]!=1) & (sus_trans[:,3]!=1) & Laser_avail & np.array([avails])).T
    
    fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]   
    fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]
    fr_sel = [f for (f, t) in zip(features_per_su, sel_avail) if t]   
    fr_nonsel = [f for (f, t) in zip(features_per_su, nonsel_avail) if t]

    coeff_all = dist_all.get_coeff(features_per_su)
    coeff_sust = dist_sust.get_coeff(fr_sust)
    coeff_trans = dist_trans.get_coeff(fr_trans)
    coeff_sel = dist_sel.get_coeff(fr_sel)
    coeff_nonsel = dist_nonsel.get_coeff(fr_nonsel)

    for r in range(repeats):
        print(r)
        (features_per_su, reg_list, _avails) = get_dataset(denovo=True, to_explore=False, random_type='one_trial')
        fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
        fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]
        fr_sel = [f for (f, t) in zip(features_per_su, sel_avail) if t]
        fr_nonsel = [f for (f, t) in zip(features_per_su, nonsel_avail) if t]

        dist_all.append_data(features_per_su, coeff_all)
        dist_sust.append_data(fr_sust, coeff_sust)
        dist_trans.append_data(fr_trans, coeff_trans)
        dist_sel.append_data(fr_sel, coeff_sel)
        dist_nonsel.append_data(fr_nonsel, coeff_nonsel)

    dist_all_list = dist_all.get_data()
    dist_sust_list = dist_sust.get_data()
    dist_trans_list = dist_trans.get_data()
    dist_sel_list = dist_sel.get_data()
    dist_nonsel_list = dist_nonsel.get_data()

    np.savez(f'D:\pixel-optogenetic\PCA\late2s\pca_dist_{repeats}.npz', dist_all=dist_all_list, dist_sust=dist_sust_list, dist_trans=dist_trans_list, dist_sel=dist_sel_list, dist_nonsel=dist_nonsel_list)
    
    with h5py.File(f'D:\pixel-optogenetic\PCA\late2s\pca_dist_{repeats}.hdf5', "w") as fw:
        fw.create_dataset('dist_sel_on', data=np.vstack(dist_sel_list[0]))   
        fw.create_dataset('dist_sel_off', data=np.vstack(dist_sel_list[1])) 
        fw.create_dataset('dist_nonsel_on', data=np.vstack(dist_nonsel_list[0]))   
        fw.create_dataset('dist_nonsel_off', data=np.vstack(dist_nonsel_list[1]))           
    
def plot_dist_one(d_one, subgroup, repeats):
    #plot time-distance 2d figure
    
    #S1-S2 distance
    
    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[0]), axis=0)
    sem = stats.sem(np.vstack(d_one[0]))
    plt.fill_between(np.arange(mm.shape[0]), mm - sem, mm + sem, color="b", alpha=0.2, ls='None')
    ax.plot(mm, '-b', lw=1)

    mm = np.mean(np.vstack(d_one[1]), axis=0)
    sem = stats.sem(np.vstack(d_one[1]))
    plt.fill_between(np.arange(mm.shape[0]), mm - sem, mm + sem, color="k", alpha=0.2, ls='None')
    ax.plot(mm, '-k', lw=1)
    ax.set_ylabel('S1-S2 PC distance')
    ax.set_xticks([12.5, 32.5, 52.5])
    ax.set_xticklabels([0, 5, 10])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [12.5, 16.5, 32.5, 34.5, 39.5, 44.5]]
    fig.savefig(f'Sample_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    #on-off distance
    # (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    # mm = np.mean(np.vstack(d_one[2]), axis=0)
    # sem = stats.sem(np.vstack(d_one[2]))
    # plt.fill_between(np.arange(mm.shape[0]), mm - sem, mm + sem, color="k", alpha=0.2, ls='None')
    # ax.plot(mm, '-k', linestyle = '--', lw=1)

    # mm = np.mean(np.vstack(d_one[3]), axis=0)
    # sem = stats.sem(np.vstack(d_one[3]))
    # plt.fill_between(np.arange(mm.shape[0]), mm - sem, mm + sem, color="k",alpha=0.2, ls='None')
    # ax.plot(mm, '-k', linestyle = ':', lw=1)
    # ax.set_ylabel('laser on-laser off PC distance')
    # ax.set_xticks([12.5, 32.5, 52.5])
    # ax.set_xticklabels([0, 5, 10])
    # ax.set_xlabel('time (s)')
    # [ax.axvline(x, color='k', ls=':') for x in [12.5, 16.5, 32.5, 34.5, 39.5, 44.5]]
    # fig.savefig(f'Laser_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')




def plot_dist(repeats):
    fstr = np.load(f'pca_dist_{repeats}.npz', allow_pickle=True)
    dist_all = fstr['dist_all']
    dist_sust = fstr['dist_sust']
    dist_trans = fstr['dist_trans']
    dist_sel = fstr['dist_sel']
    dist_nonsel = fstr['dist_nonsel']
    
    plot_dist_one(dist_all, 'all', repeats)
    plot_dist_one(dist_sust, 'sust', repeats)
    plot_dist_one(dist_trans, 'trans', repeats)
    plot_dist_one(dist_sel, 'sel', repeats)
    plot_dist_one(dist_nonsel, 'nosel', repeats)


if __name__ == "__main__":
    # rcParams['pdf.fonttype'] = 42
    # rcParams['ps.fonttype'] = 42
    # rcParams['font.family'] = 'sans-serif'
    # rcParams['font.sans-serif'] = ['Arial']    
    plot3d_all(plot_delay=6, repeats=100)
    dist_all(repeats=100)
    # plot_dist(repeats=500)
    # delay = 3
    # (features_per_su, reg_list) = get_dataset(True, rnd_half=True)
    # trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    # # list(trans6.keys())
    # sust = trans_fstr["sust"]
    # trans = trans_fstr["transient"]
