#!/usr/bin/env python
import os
import sys
import glob
import multiprocessing as mp
import numpy as np
from sklearn.cluster import AgglomerativeClustering

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='''TrRefine: Refinement of trRosetta outputs''')

    # Input npz file (outputs of msa-net, tbm-net)
    parser.add_argument("-msa_npz", required=True, \
                         help="output npz file of trRosetta msa-net")
    parser.add_argument("-tbm_npz", required=True, \
                         help="output npz file of trRosetta tbm-net")
    parser.add_argument("-pdb_dir_s", required=True, nargs='+', \
                         help="path to predicted pdb files by trRosetta")
    parser.add_argument("-a3m_fn", required=True, \
                         help="MSA file in a3m format")
    parser.add_argument("-hhr_fn", required=True, \
                         help="MSA file in hhr format")
    parser.add_argument('-n_core', required=True, type=int,\
                         help='The number of cores can ben used')
    parser.add_argument("--rescore", action='store_true', default=False)
    args = parser.parse_args()

    return args

def rescore(pose):
    score_fxn = pyrosetta.create_score_function('ref2015_cart')
    score = score_fxn.score(pose)
    return score

def calc_lddt_dist(args):
    i, j, pose_s = args
    pose_i = pose_s[i]
    pose_j = pose_s[j]
    #
    lddt_1 = float(os.popen("/home/minkbaek/bin/lddt -c %s %s | grep Glob"%(pose_i, pose_j)).readlines()[-1].split()[-1])
    lddt_2 = float(os.popen("/home/minkbaek/bin/lddt -c %s %s | grep Glob"%(pose_j, pose_i)).readlines()[-1].split()[-1])
    lddt = (lddt_1 + lddt_2) / 2.0
    return 1 - lddt

def pick_rep(msa_npz, tbm_npz, pdb_dir_s, n_core, n_clust=10, rescore=False, cutoff=0.8):
    # pick 10 lowest E conformation from each method
    pdb_fn_s = list()
    score_s = list()
    pose_s = list()
    for pdb_dir in pdb_dir_s:
        fn_s = glob.glob("%s/model*.pdb"%(pdb_dir))
        fn_s.sort()
        #
        if not rescore:
            scores = list()
            for fn in fn_s:
                sc = float(os.popen("grep ^pose %s"%fn).readlines()[-1].split()[-1])
                scores.append(sc)
        else:
            import pyrosetta
            pyrosetta.init('-mute all')

            poses = list()
            for pdb_fn in fn_s:
                pose = pyrosetta.pose_from_file(pdb_fn)
                poses.append(pose.clone())

            # Setup multiprocessor
            print ("setup multiprocessor")
            n_core_pool = min(n_core, len(fn_s))
            pool = mp.Pool(n_core_pool)
            #
            # rescore all the inputs using ref2015_cart & filter out highE conf
            print ("rescoring")
            scores = pool.map(rescore, poses)
            pool.close()
            pool.join()
        scores = np.array(scores)
        idx = np.argsort(scores)
        fn_s = np.array(fn_s)
        pdb_fn_s.append(fn_s[idx[:10]])
        score_s.append(scores[idx[:10]])
        sys.stdout.write("INFO: Pick 10 lowest E conformations from %s, Emin=%.3f / Emax=%.3f\n"%(pdb_dir, score_s[-1][0], score_s[-1][-1]))

    pdb_fn_s = np.concatenate(pdb_fn_s)
    score_s = np.concatenate(score_s)
    #
    Emin = np.min(score_s)
    Ecut = cutoff * Emin
    idx_s = np.where(score_s < Ecut)[0]
    n_clust = min(len(idx_s), n_clust)
    filtered = list()
    for idx in idx_s:
        filtered.append(pdb_fn_s[idx])
    score_s = score_s[idx_s]
    model_s = np.array(pdb_fn_s)[idx_s]
    n_str = len(filtered)
    #
    sys.stdout.write("INFO: After filtering based on E, Emin=%.3f Ecut=%.3f, %d / %d\n"%(Emin, Ecut, len(filtered), len(pdb_fn_s)))
    del pdb_fn_s

    args = list()
    for i in range(n_str-1):
        for j in range(i+1, n_str):
            args.append((i,j,filtered))

    n_core_pool = min(n_core, len(args))
    pool = mp.Pool(n_core_pool)
    raw_dist = pool.map(calc_lddt_dist, args)
    pool.close()
    pool.join()
    dist = np.zeros((n_str, n_str), dtype=np.float)
    idx = np.triu_indices(n_str, k=1)
    dist[idx] = raw_dist
    dist = dist + dist.T
    #
    cluster = AgglomerativeClustering(n_clusters=n_clust, affinity='precomputed', linkage='single').fit(dist)
    #
    unique_labels = np.unique(cluster.labels_)
    rep_s = list()
    for label in unique_labels:
        idx = np.where(cluster.labels_==label)[0]
        #
        Emin_idx = np.argmin(score_s[idx])
        model = model_s[idx][Emin_idx]
        rep_s.append(os.path.abspath(model))
    #
    if not os.path.exists("rep_s"):
        os.mkdir("rep_s")
    os.chdir('rep_s')
    input_s = list()
    distogram_s = list()
    for i_rep, rep in enumerate(rep_s):
        os.system("ln -sf %s rep_%d.pdb"%(rep, i_rep))
        if "pdb-msa" in rep:
            distogram_s.append(msa_npz)
        else:
            distogram_s.append(tbm_npz)
        input_s.append(os.path.abspath("rep_%d.pdb"%i_rep))
    #
    with open('inpdb_dan.list', 'wt') as fp:
        fp.write("\n".join(input_s))
        fp.write("\n")
    with open('distogram_dan.list', 'wt') as fp:
        fp.write("\n".join(distogram_s))
        fp.write("\n")

    os.chdir('..')
    return rep_s

def main():
    args = get_args()

    # Pick 10 representative structures from the given pdbs
    if not os.path.exists("rep_s/distogram_dan.list"):
        rep_s = pick_rep(args.msa_npz, args.tbm_npz, args.pdb_dir_s, args.n_core, rescore=args.rescore)
    else:
        rep_s = [line.strip() for line in open("rep_s/inpdb_dan.list")]

    # run DAN-msa
    script_dir = os.path.dirname(__file__)
    acc_fn_s = glob.glob("rep_s/*_acc.npz")
    if len(acc_fn_s) < len(rep_s):
        os.system("python -u -W ignore %s/DAN-msa/ErrorPredictorMSA.py -p %d rep_s/distogram_dan.list rep_s/inpdb_dan.list rep_s"%(script_dir, args.n_core))

    # run SStor
    SS_fn_s = glob.glob("rep_s/*_SS.npz")
    if len(SS_fn_s) < len(rep_s):
        os.system("python -u -W ignore %s/SStor_pred/main_multi.py -a3m_fn %s -hhr_fn %s -pdb_fn_s rep_s/rep_?.pdb"%(script_dir, args.a3m_fn, args.hhr_fn))
    #
    tor_s = list()
    for i in range(len(rep_s)):
        tor_s.append(np.load("rep_s/rep_%d_SS.npz"%i)['tor'].astype(np.float32)+1e-9)
    tor_s = np.stack(tor_s, axis=0)
    tor_s = np.mean(tor_s, axis=0)
    np.savez_compressed("rep_s/BBtor.npz", phi=tor_s[:,:36], psi=tor_s[:,36:72], omega=tor_s[:,72:])
        
    
    # run trRefine
    if not os.path.exists("t000_.trRefine.npz"):
        os.system("python -u -W ignore %s/trRefine/main.py -a3m_fn %s -pdb_fn_s rep_s/rep_?.pdb -acc_fn_s rep_s/rep_?_acc.npz -SS_fn_s rep_s/rep_?_SS.npz -npz_fn_s %s %s -out_fn t000_.trRefine.npz"%(script_dir, args.a3m_fn, args.msa_npz, args.tbm_npz))

main()
