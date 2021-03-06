#!/usr/bin/python

import warnings, logging, os, sys
warnings.filterwarnings('ignore',category=FutureWarning)
logging.disable(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import json
import tensorflow as tf
import numpy as np
import string

if len(sys.argv) < 3:
    print ("python %s [MSA file (a3m or FASTA format)] [output npz file]"%__file__)
    sys.exit()

msa_file = sys.argv[1]
npz_file = sys.argv[2]

MDIR     = '%s/../weights/trrosetta_homo'%(os.path.dirname(os.path.realpath(__file__)))

n2d_layers   = 36
n2d_filters  = 64
window2d     = 3
wmin         = 0.8
ns           = 21

# read A3M and convert letters into
# integers in the 0..20 range
def parse_a3m(filename):
    seqs = []
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    # read file line by line
    for line in open(filename,"r"):
        # skip labels
        if len(line.strip()) == 0:
            continue
        if line[0] != '>':
            # remove lowercase letters and right whitespaces
            seqs.append(line.rstrip().translate(table))

    # convert letters into numbers
    alphabet = np.array(list("ARNDCQEGHILKMFPSTWYV-"), dtype='|S1').view(np.uint8)
    msa = np.array([list(s) for s in seqs], dtype='|S1').view(np.uint8)
    for i in range(alphabet.shape[0]):
        msa[msa == alphabet[i]] = i

    # treat all unknown characters as gaps
    msa[msa > 20] = 20
    return msa

# 1-hot MSA to PSSM
def msa2pssm(msa1hot, w):
    beff = tf.reduce_sum(w)
    f_i = tf.reduce_sum(w[:,None,None]*msa1hot, axis=0) / beff + 1e-9
    h_i = tf.reduce_sum( -f_i * tf.math.log(f_i), axis=1)
    return tf.concat([f_i, h_i[:,None]], axis=1)

# reweight MSA based on cutoff
def reweight(msa1hot, cutoff):
    with tf.name_scope('reweight'):
        id_min = tf.cast(tf.shape(msa1hot)[1], tf.float32) * cutoff
        id_mtx = tf.tensordot(msa1hot, msa1hot, [[1,2], [1,2]])
        id_mask = id_mtx > id_min
        w = 1.0/tf.reduce_sum(tf.cast(id_mask, dtype=tf.float32),-1)
    return w

# shrunk covariance inversion
def fast_dca(msa1hot, weights, penalty = 4.5):
    nr = tf.shape(msa1hot)[0]
    nc = tf.shape(msa1hot)[1]
    ns = tf.shape(msa1hot)[2]
    with tf.name_scope('covariance'):
        x = tf.reshape(msa1hot, (nr, nc * ns))
        num_points = tf.reduce_sum(weights) - tf.sqrt(tf.reduce_mean(weights))
        mean = tf.reduce_sum(x * weights[:,None], axis=0, keepdims=True) / num_points
        x = (x - mean) * tf.sqrt(weights[:,None])
        cov = tf.matmul(tf.transpose(x), x)/num_points
    with tf.name_scope('inv_convariance'):
        cov_reg = cov + tf.eye(nc * ns) * penalty / tf.sqrt(tf.reduce_sum(weights))
        inv_cov = tf.linalg.inv(cov_reg)
        x1 = tf.reshape(inv_cov,(nc, ns, nc, ns))
        x2 = tf.transpose(x1, [0,2,1,3])
        features = tf.reshape(x2, (nc, nc, ns * ns))
        x3 = tf.sqrt(tf.reduce_sum(tf.square(x1[:,:-1,:,:-1]),(1,3))) * (1-tf.eye(nc))
        apc = tf.reduce_sum(x3,0,keepdims=True) * tf.reduce_sum(x3,1,keepdims=True) / tf.reduce_sum(x3)
        contacts = (x3 - apc) * (1-tf.eye(nc))
    return tf.concat([features, contacts[:,:,None]], axis=2)

a3m = parse_a3m(msa_file)

contacts = {'pd':[], 'ph':[], 'po':[], 'pt':[], 'pp':[]}

#
# network
#
config = tf.ConfigProto(
    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.9)
)
activation = tf.nn.elu
conv1d = tf.layers.conv1d
conv2d = tf.layers.conv2d
with tf.Graph().as_default():

    with tf.name_scope('input'):
        ncol = tf.placeholder(dtype=tf.int32, shape=())
        nrow = tf.placeholder(dtype=tf.int32, shape=())
        msa = tf.placeholder(dtype=tf.uint8, shape=(None,None))
        is_train = tf.placeholder(tf.bool, name='is_train')

    #
    # collect features
    #
    msa1hot  = tf.one_hot(msa, ns, dtype=tf.float32)
    w = reweight(msa1hot, wmin)

    # 1D features
    f1d_seq = msa1hot[0,:,:20]
    f1d_pssm = msa2pssm(msa1hot, w)

    f1d = tf.concat(values=[f1d_seq, f1d_pssm], axis=1)
    f1d = tf.expand_dims(f1d, axis=0)
    f1d = tf.reshape(f1d, [1,ncol,42])

    # 2D features
    f2d_dca = tf.cond(nrow>1, lambda: fast_dca(msa1hot, w), lambda: tf.zeros([ncol,ncol,442], tf.float32))
    f2d_dca = tf.expand_dims(f2d_dca, axis=0)

    f2d = tf.concat([tf.tile(f1d[:,:,None,:], [1,1,ncol,1]),
                    tf.tile(f1d[:,None,:,:], [1,ncol,1,1]),
                    f2d_dca], axis=-1)
    f2d = tf.reshape(f2d, [1,ncol,ncol,442+2*42])


    #
    # 2D network
    #
    layers2d = [f2d]
    layers2d.append(conv2d(layers2d[-1], n2d_filters, 1, padding='SAME'))
    layers2d.append(tf.contrib.layers.instance_norm(layers2d[-1]))
    layers2d.append(activation(layers2d[-1]))

    # stack of residual blocks with dilations
    dilation = 1
    for _ in range(n2d_layers):
        layers2d.append(conv2d(layers2d[-1], n2d_filters, window2d, padding='SAME', dilation_rate=dilation))
        layers2d.append(tf.contrib.layers.instance_norm(layers2d[-1]))
        layers2d.append(activation(layers2d[-1]))
        layers2d.append(tf.keras.layers.Dropout(rate=0.15)(layers2d[-1], training=is_train))
        layers2d.append(conv2d(layers2d[-1], n2d_filters, window2d, padding='SAME', dilation_rate=dilation))
        layers2d.append(tf.contrib.layers.instance_norm(layers2d[-1]))
        layers2d.append(activation(layers2d[-1] + layers2d[-7]))
        dilation *= 2
        if dilation > 16:
            dilation = 1

    # anglegrams for theta
    logits_theta = conv2d(layers2d[-1], 25, 1, padding='SAME')
    prob_theta = tf.nn.softmax(logits_theta)

    # anglegrams for phi
    logits_phi = conv2d(layers2d[-1], 13, 1, padding='SAME')
    prob_phi = tf.nn.softmax(logits_phi)

    # symmetrize
    layers2d.append(0.5 * (layers2d[-1] + tf.transpose(layers2d[-1], perm=[0,2,1,3])))

    # distograms
    logits_dist = conv2d(layers2d[-1], 37, 1, padding='SAME')
    prob_dist = tf.nn.softmax(logits_dist)

    # beta-strand pairings (not used)
    logits_bb = conv2d(layers2d[-1], 3, 1, padding='SAME')
    prob_bb = tf.nn.softmax(logits_bb)

    # anglegrams for omega
    logits_omega = conv2d(layers2d[-1], 25, 1, padding='SAME')
    prob_omega = tf.nn.softmax(logits_omega)

    # hoo-distograms
    logits_hdist = conv2d(layers2d[-1], 37, 1, padding='SAME')
    prob_hdist = tf.nn.softmax(logits_hdist)

    saver = tf.train.Saver()

#    for ckpt in ['default_xaa2.meta', 'default_xaa2.meta', 'default_xaa2.meta']:
    for filename in os.listdir(MDIR):
        if not filename.endswith(".index"):
            continue
        ckpt = MDIR+"/"+os.path.splitext(filename)[0]
        print("CHECKPOINT: " + ckpt + "\n\n\n\n")
        with tf.Session(config=config) as sess:
            saver.restore(sess, ckpt)
            pd, ph, pt, pp, po = sess.run([prob_dist, prob_hdist, prob_theta, prob_phi, prob_omega],
                                       feed_dict = {msa : a3m, ncol : a3m.shape[1], nrow : a3m.shape[0], is_train : 0})
            contacts['pd'].append(pd[0])
            contacts['ph'].append(ph[0])
            contacts['pt'].append(pt[0])
            contacts['po'].append(po[0])
            contacts['pp'].append(pp[0])
            print(ckpt, '- done')

# average over all network params
contacts['pd'] = np.mean(contacts['pd'], axis=0)
contacts['ph'] = np.mean(contacts['ph'], axis=0)[:,:,1:21].sum(axis=-1)
contacts['pt'] = np.mean(contacts['pt'], axis=0)
contacts['po'] = np.mean(contacts['po'], axis=0)
contacts['pp'] = np.mean(contacts['pp'], axis=0)

# save distograms & anglegrams
np.savez_compressed(npz_file, dist=contacts['pd'], hcont=contacts['ph'], omega=contacts['po'], theta=contacts['pt'], phi=contacts['pp'])
