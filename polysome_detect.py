import starfile
import numpy as np
import pandas as pd
from disjoint_set import DisjointSet
from tqdm import tqdm
from matplotlib import pyplot as plt
from fromEuler_RELION import fromEuler
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--threshold", required=True, type=float, help="threshold for ribosomes considered to be in polysomes, in pixels")

    args = parser.parse_args()

    pixel_threshold = args.threshold

    # read star file
    df = starfile.read("refine080_data.star")["particles"]

    # set values
    boxsize=np.array([256, 256, 256])
    # TODO: WHAT IS PIXEL SIZE?
    pixelsize=1.7005
    localxyz_mRNAentry= np.array([128.721, 241.009, 236.059])/pixelsize
    offset_mRNAentry = localxyz_mRNAentry-(boxsize+1)/2
    localxyz_mRNAexit =np.array([186.995, 248.341, 300])/pixelsize
    offset_mRNAexit = localxyz_mRNAexit-(boxsize+1)/2

    n = df.shape[0]
    coords = df.iloc[:,0:3].to_numpy()
    origins = df.iloc[:,18:21].to_numpy()
    angles = df.iloc[:,3:6].to_numpy()

    R = np.zeros([n,3,3])
    for i in range(n):
        rot = angles[i,0]
        tilt = angles[i,1]
        psi = angles[i,2]

        R[i] = fromEuler(rot, tilt, psi)

    # print(R)
    # print(R.shape)

    # TODO: FIGURE OUT UNITS
    entry_exit = np.zeros([n,6])
    entry_exit[:,0:3] = coords - origins + np.dot(R, offset_mRNAentry)
    entry_exit[:,3:6] = coords - origins + np.dot(R, offset_mRNAexit)

    np.savetxt("entry_exit.csv", entry_exit, delimiter=',')

    # detect polysomes
    tomoNames = df.loc[:, "rlnMicrographName"]
    uniqueTomoNames = pd.unique(tomoNames)
    tomoNum = uniqueTomoNames.shape[0]

    counter = 0
    angle_differences = []
    for t in tqdm(range(tomoNum)):
        index = df.index[df["rlnMicrographName"] == uniqueTomoNames[t]][0]

        particles = df.loc[df["rlnMicrographName"] == uniqueTomoNames[t]]
        pn = particles.shape[0]

        particles_np = particles.to_numpy()
        # polysomes = np.zeros([pn, 2])

        # TODO: find polysomes with rank information
        # polyid 0: monosome, >0: belongs to polysome
        # polyid = 0 

        ds = DisjointSet()

        distances = np.zeros([pn,pn])
        for p1 in range(pn):
            ds.find(p1)
            for p2 in range(pn):
                if p1 != p2:
                    ds.find(p2)
                    distance = np.linalg.norm(entry_exit[p1,3:6] - entry_exit[p2,0:3])
                    distances[p1,p2] = distance
                    # threshold is 41 pixels = 7 nm
                    # try 410?
                    if distance < pixel_threshold:
                        ds.union(p1, p2)

        np.savetxt("distances.csv", distances, delimiter=',')
        polysomes = list(ds.itersets())
        dimers = [x for x in polysomes if len(x) == 2]

        for dimer in dimers:
            dimer = list(dimer)
            r1 = index + dimer[0]
            r2 = index + dimer[1]

            r1_psi = angles[r1,2]
            r1_phi = angles[r1,1]
            r2_psi = angles[r2,2]
            r2_phi = angles[r2,1]

            angle_differences.append((r2_phi - r1_phi, r2_psi - r1_psi))

        # print(dimers)
        # print(pn)

        counter += 1
        
    angle_differences_np = np.array(angle_differences)
    plt.scatter(angle_differences_np[:,0],angle_differences_np[:,1])
    plt.savefig("angle_differences_{}.png".format(pixel_threshold))
    # print(angle_differences)