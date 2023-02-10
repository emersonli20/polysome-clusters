import numpy as np

def fromEuler(rot, tilt, psi):
    ca = np.cos(rot);
    cb = np.cos(tilt);
    cg = np.cos(psi);
    sa = np.sin(rot);
    sb = np.sin(tilt);
    sg = np.sin(psi);
    cc = cb * ca;
    cs = cb * sa;
    sc = sb * ca;
    ss = sb * sa;

    A = np.zeros([3,3])
    A[0,0] = cg * cc - sg * sa;
    A[0,1] = cg * cs + sg * ca;
    A[0,2] = -cg * sb;
    A[1,0] = -sg * cc - cg * sa;
    A[1,1] = -sg * cs + cg * ca;
    A[1,2] = sg * sb;
    A[2,0] = sc;
    A[2,1] = ss;
    A[2,2] = cb;

    return A