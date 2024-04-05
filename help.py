import numpy as np
from scipy.integrate import trapezoid

def SQ(file, size = -1):

    qx = np.loadtxt(file, skiprows=1, max_rows=1)
    Stotal = np.loadtxt(file, skiprows=3)
    maxQ = qx[-1]
    qn = len(qx)
    n = np.shape(Stotal)[0] // qn
    Stotal = np.reshape(Stotal, (n, qn, qn))
    
    if size == -1:
        size = len(qx)//2
    else: 
        size = int(size)
    finalS = np.zeros(size)
    
    for S in Stotal:
        S[(qn - 1)//2, (qn -1)//2] = np.min(S)
        index = (np.sqrt(qx[:, np.newaxis]**2 + qx**2) / maxQ * size).astype(int)
        index[index >= size] = size - 1
        count = np.bincount(index.flatten(), minlength=size)
        s = np.bincount(index.flatten(), weights=S.flatten())  
    
        finalS += s/count
    
    qq = np.sort(np.abs(qx))
    q = np.linspace(qq[1], qq[-1], size)
    return q[1:], finalS[1:]/n 

def integrate(t, F, w):
    return trapezoid(F*np.cos(w*t),t)

def FQT(file, size = -1, N = 5000):
    qx = np.loadtxt(file, skiprows=1, max_rows=1)
    t = np.loadtxt(file, skiprows=2, max_rows=1)
    Stotal = np.loadtxt(file, skiprows=3, dtype = np.complex_)
    maxQ = qx[-1]
    qn = len(qx)
    n = np.shape(Stotal)[0] // qn
    if size == -1:
        size = len(qx)//2
    else: 
        size = int(size)
    FQT = np.reshape(Stotal, (n, qn, qn))
    FQTav = np.zeros((n, size))
    index = (np.sqrt(qx[:, np.newaxis]**2 + qx**2) / maxQ * size).astype(int)
    index[index >= size] = size - 1
    count = np.bincount(index.flatten(), minlength=size)
    for i in range(n):
        FQTav[i] = np.bincount(index.flatten(), weights=FQT[i].real.flatten(), minlength=size)/count
    return t, qx, FQTav[:, 1:]/FQTav[0, 1:]