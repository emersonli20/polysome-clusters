import starfile
import numpy as np
import pandas as pd

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

# TODO: use the values above to calculate the entry and exit sites of every particle

# detect polysomes
tomoNames = df.loc[:, "rlnMicrographName"]
uniqueTomoNames = pd.unique(tomoNames)
tomoNum = uniqueTomoNames.shape[0]

for t in range(tomoNum):
    particles = df.loc[df["rlnMicrographName"] == uniqueTomoNames[t]]
    pn = particles.shape[0]

    for p in range(pn):
        for p2 in range(pn):
            if p != p2:
                # TODO: calculate distance from exit to entry site
                distance = np.linalg.norm()
    print(pn)
    break
    