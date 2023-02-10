import starfile
import numpy as np
import pandas as pd
from fromEuler_RELION import fromEuler

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

print(R)
print(R.shape)

# TODO: FIGURE OUT UNITS
entry_exit = np.zeros([n,6])
entry_exit[:,0:3] = coords - origins + np.dot(R, offset_mRNAentry)
entry_exit[:,3:6] = coords - origins + np.dot(R, offset_mRNAexit)

# detect polysomes
tomoNames = df.loc[:, "rlnMicrographName"]
uniqueTomoNames = pd.unique(tomoNames)
tomoNum = uniqueTomoNames.shape[0]

for t in range(tomoNum):
    particles = df.loc[df["rlnMicrographName"] == uniqueTomoNames[t]]
    pn = particles.shape[0]

    particles_np = particles.to_numpy()
    polysomes = np.zeros([pn, 2])
    # polyid 0: monosome, >0: belongs to polysome
    polyid = 0 

    for p in range(pn):
        for p2 in range(pn):
            if p != p2:
                # TODO: calculate distance from exit to entry site
                distance = np.linalg.norm(entry_exit[:,3:6] - entry_exit[:,0,3])
                # threshold is 41 pixels = 7 nm
                if distance < 41:
                    # case 1: p and p2 are both in polysome already
                    if polysomes[p,0] > 0 and polysomes[p2,0] > 0:
                        # TODO: CONVERT NEXT 2 LINES TO PYTHON
                        # motl(7,motl(6,:)==motl(6,p2)) = motl(7,motl(6,:)==motl(6,p2)) + motl(7,p);
                        # motl(6,motl(6,:)==motl(6,p2)) = motl(6,p);  
                    # case 2: only p2 is in polysome already
                    if polysomes[p,0] == 0 and polysomes[p2,0] > 0:
                        polysomes[p,0] = polysomes[p2,0]
                        # TODO: CONVERT NEXT LINE TO PYTHON
                        # motl(7,motl(6,:)==motl(6,p2)) = motl(7,motl(6,:)==motl(6,p2)) +1;                         
                    # case 3: only p1 is in polysome already
                    if polysomes[p,0] > 0 and polysomes[p2,0] == 0:
                        polysomes[p2,0] = polysomes[p,0]
                        polysomes[p2,1] = polysomes[p2,0] + 1
                    # case 4: neither is in polysome already
                    if polysomes[p,0] == 0 and polysomes[p2,0] == 0:
                        polyid = polyid + 1
                        polysomes[p,0] = polyid
                        polysomes[p2,0] = polyid
                        polysomes[p,1] = 1
                        polysomes[p2,1] = 2
                        
    print(pn)
    break
    