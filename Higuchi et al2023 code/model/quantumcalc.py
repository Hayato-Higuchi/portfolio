import numpy as np
from qiskit import *
import scipy.linalg as LA
from itertools import product
from qiskit.quantum_info.operators import Operator
from qiskit.circuit.library import MCXGate

def generate_Bmat(A):
    matsize = np.shape(A)[0]
    ar = np.identity(matsize) - A.dot(A)
    B1 = A+1.0j*LA.sqrtm(ar)
    B2 = A-1.0j*LA.sqrtm(ar)
    
    B =  np.zeros((2*matsize,2*matsize),dtype=np.complex128)
    
    for i in range(2*matsize):
        if i//matsize == 0:
            B[i,i] = B1[i,i]
        else:
            B[i,i] = B2[matsize-i,matsize-i]
    
    return B

def generate_CB_gate(num_c,Bmat,state=''):
    
    if len(state) == 0:
        state = '1'*num_c
        
    matrixsize = np.shape(Bmat)[0]
    
    num_t = int(np.log2(matrixsize))
    N = num_c + num_t
    
    ID = np.identity(2**num_t)
    
    CB = np.zeros((2**N,2**N),dtype=np.complex128)

    for i in range(2**num_c):
        ibin = bin(i)[2:].zfill(num_c)
        projection = np.zeros((2**num_c,2**num_c),dtype=np.complex128)
        projection[i,i] = 1.0+0.0j
        if ibin == state:
            CB += np.kron(Bmat,projection)
        else:
            CB += np.kron(ID,projection)
    
    ctr = QuantumRegister(num_c,name='ctr')
    trg = QuantumRegister(num_t,name='trg')
    QC = QuantumCircuit(ctr,trg,name=f'CB_{state}')
    CBO = Operator(CB)
    CBO.name = 'controlled-B'
    QC.append(CBO,range(N))
    
    return QC

def Nloops(X,operation_name='ret'):
    let = f'{operation_name} = product('
    for x in X:
        let += f'range({x}),'
    let = let[:-1] + ')'

    return let

def circuitize(operator,operator_name=''):
    site = int(np.log2(np.shape(operator)[0]))
    
    QC = QuantumCircuit(site,name=operator_name)
    Op = Operator(operator)
    QC.append(Op,range(site))
    
    return QC

def check_incre_dcre(I,D):
    Ip = Operator(I).data
    Dp = Operator(D).data
    N = np.shape(Ip)[0]
    
    return np.allclose(np.dot(Ip,Dp),np.identity(N))

def MCX_with_state(num_c,state=''):
    if not len(state):
        state = '1'*num_c
    if len(state) != num_c:
        print('STATE DESIGNATION ERROR')
        return None
    QC = QuantumCircuit(num_c+1,name=f'MCX_{state}')
    for i in range(num_c):
        if state[i] == '0':
            QC.x(i)
    QC.append(MCXGate(num_c),range(num_c+1))
    for i in range(num_c):
        if state[i] == '0':
            QC.x(i)
    return QC

def normalized(F,fac=[]):
    f_norm = np.sqrt(sum([ f_**2 for f_ in F]))
    fac.append(f_norm)
    return np.array([f/f_norm for f in F])

def flatten_and_normalized(F):
    F = F.flatten()
    f_norm = np.sqrt(sum([ f_**2 for f_ in F]))
    return np.array([f/f_norm for f in F])

def incre(n):
    QC = QuantumCircuit(n)
    for i in range(n-1):
        QC.append(MCXGate(n-1-i),range(n-i))
    QC.x(0)
    return QC

def decre(n):
    QC = QuantumCircuit(n)
    for i in range(n-1):
        for ii in range(n-1-i):
            QC.x(n-2-i-ii)
        QC.append(MCXGate(n-1-i),range(n-i))
        for ii in range(n-1-i):
            QC.x(n-2-i-ii)
    QC.x(0)
    return QC

def control_incre(num_ctrl,num_target,state):
    QC = QuantumCircuit(num_ctrl+num_target,name=f'Ic({state[::-1]})')
    Inc = incre(num_target)
    Ic = Inc.control(num_ctrl)
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    QC.append(Ic,range(num_ctrl+num_target))
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    return QC

def control_decre(num_ctrl,num_target,state):
    QC = QuantumCircuit(num_ctrl+num_target,name=f'Dc({state[::-1]})')
    Dec = decre(num_target)
    Dc = Dec.control(num_ctrl)
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    QC.append(Dc,range(num_ctrl+num_target))
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    return QC

def control_hadamard(num_ctrl,state):
    QC = QuantumCircuit(num_ctrl+1,name=f'Hc')#({state})')
    qch = QuantumCircuit(1)
    qch.h(0)
    Hc = qch.control(num_ctrl)
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    QC.append(Hc,range(num_ctrl+1))
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    return QC

def applying_B(X, QC):
    state = Operator(QC).data[:,0]
    for i in range(np.shape(state)[0]):
        state[i] *= X[i]
    state = normalized(state)
    qc = QuantumCircuit(QC.width())
    qc.initialize(state,range(QC.width()))
    return qc

def apply_Bope(coefs,qc,sites,anc,fac):
    N = len(sites)
    A = normalized(coefs,fac)
    B = np.zeros((2**(N+1),2**(N+1)),dtype=np.float64())
    matrixsize = 2**(N+1)
    M = 2**(N)
    for i in range(matrixsize):
        for j in range(matrixsize):
            if i == j and i < M:
                B[i,i] = A[i]
            elif i == j and i >= M:
                B[i,i] = - A[i-M]
            elif i == j+ M:
                B[i,j] = np.sqrt( 1 - A[j]**2)
            elif j == i+ M:
                B[i,j] = np.sqrt( 1 - A[i]**2)
    B_ope = Operator(B)
    qc.append(B_ope,sites+[anc])
    return qc

def control_inversion(num_ctrl,state):
    QC = QuantumCircuit(num_ctrl+2,name=f'Inv_c({state})')
    qccz = QuantumCircuit(2)
    qccz.cz(0,1)
    Inv_c = qccz.control(num_ctrl)
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    QC.append(Inv_c,range(num_ctrl+2))
    for i in range(num_ctrl):
        if state[i] == '0':
            QC.x(i)
    return QC

def get_Maxwell_res(state,species,direction,qubit_mesh):
    if species == 'E':
        spe_ind = '0'
    elif species == 'B':
        spe_ind = '1'
    else:
        print('wrong species')
        
    if direction == 'x':
        dir_ind = '00'
    elif direction == 'y':
        dir_ind = '01'
    elif direction == 'z':
        dir_ind = '10'
    else:
        print('wrong direction')
    # ancila/sub/direction/species/timeの順

    ret = np.zeros(2**(3*qubit_mesh))
    for i in range(len(state)):
        ibin = bin(i)[2:].zfill(int(np.log2(len(state))))
        ibin_sub = ibin[:8]
        ibin_phys = ibin[8:]
        if ibin_sub == '0' + '000' + dir_ind + spe_ind + '0':
            ret[int(ibin_phys,2)] = state[i]
            
    return ret
