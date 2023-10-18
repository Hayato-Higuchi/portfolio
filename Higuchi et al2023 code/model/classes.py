import numpy as np
from qiskit import *

class Mesh():
    def __init__(self,x,y,z,vx,vy,vz):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
    
    def data(self):
        return [self.x,self.y,self.z,self.vx,self.vy,self.vz]

class Mesh_m():
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    def data(self):
        return [self.x,self.y,self.z]

class Geometry():
    def __init__(self,x,y,z,vx,vy,vz):
        self.lx = 2**x
        self.ly = 2**y
        self.lz = 2**z
        self.lvx = 2**vx
        self.lvy = 2**vy
        self.lvz = 2**vz
        self.volume = 2**x * 2**y * 2**z * 2**vx * 2**vy * 2**vz
    
    def length(self):
        return [self.lx,self.ly,self.lz,self.lvx,self.lvy,self.lvz]
    
    def grids(self,mins,maxs):
        ret = []
        Gls = self.length()
        for i in range(6):
            mi = mins[i]
            ma = maxs[i]
            gl = Gls[i]
            ret.append(np.linspace(mi,ma,gl) )
        return ret
    
class Geometry_m():
    def __init__(self,x,y,z):
        self.lx = 2**x
        self.ly = 2**y
        self.lz = 2**z
        self.volume = 2**x * 2**y * 2**z
    
    def length(self):
        return [self.lx,self.ly,self.lz]
    
    def grids(self,mins,maxs):
        ret = []
        Gls = self.length()
        for i in range(3):
            mi = mins[i]
            ma = maxs[i]
            gl = Gls[i]
            ret.append(np.linspace(mi,ma,gl) )
        return ret

class Number_of_qubits():
    def __init__(self,n_phys,n_order,n_species,n_subnode,n_ancila):
        self.ancila = n_ancila
        self.sub = n_subnode
        self.order = n_order
        self.species = n_species
        self.phys = n_phys
        self.total = n_ancila + n_phys + n_subnode + n_order + n_species
        self.control = n_subnode + n_ancila
        self.data = [n_phys,n_order,n_species,n_subnode,n_ancila]

class Number_of_qubits_m():
    def __init__(self,n_phys,n_time,n_species,n_dir,n_subnode,n_ancila):
        self.ancila = n_ancila
        self.sub = n_subnode
        self.dir = n_dir
        self.species = n_species
        self.time = n_time
        self.phys = n_phys
        self.total = n_ancila + n_phys + n_time + n_dir + n_species + n_subnode 
        self.control = n_time + n_dir + n_species + n_subnode 
        self.data = [n_phys,n_time,n_species,n_dir,n_subnode,n_ancila]

class Number_of_qubits_sampling():
    def __init__(self,n_phys,n_order,n_species,n_subnode,n_ancila,n_cla):
        self.classic = n_cla
        self.ancila = n_ancila
        self.sub = n_subnode
        self.order = n_order
        self.species = n_species
        self.phys = n_phys
        self.total = n_ancila + n_phys + n_subnode + n_order + n_species
        self.control = n_subnode + n_ancila
        self.data = [n_phys,n_order,n_species,n_subnode,n_ancila,n_cla]

class Calculation_circuit():
    def __init__(self,Ndata):
        self.n_phys = Ndata[0]
        self.n_order = Ndata[1]
        self.n_species = Ndata[2]
        self.n_subnode = Ndata[3]
        self.n_ancila = Ndata[4]
    
    def generate(self,registers):

        Rphys = QuantumRegister(self.n_phys,name='phys')
        Rorder = QuantumRegister(self.n_order,name='order')
        Rspecies = QuantumRegister(self.n_species,name='species')
        Rsub = QuantumRegister(self.n_subnode,name='sub')
        Ra = QuantumRegister(self.n_ancila,name='a')
        
        QC = QuantumCircuit(Rphys,Rorder,Rspecies,Rsub,Ra)
        
        registers.set_phys(Rphys[::-1])
        registers.set_order(Rorder[::-1])
        registers.set_species(Rspecies[::-1])
        registers.set_sub(Rsub[::-1])
        registers.set_a(Ra[::-1])
        return QC

class Calculation_circuit_m():
    def __init__(self,Ndata):
        self.n_phys = Ndata[0]
        self.n_time = Ndata[1]
        self.n_species = Ndata[2]
        self.n_dir = Ndata[3]
        self.n_subnode = Ndata[4]
        self.n_ancila = Ndata[5]
    
    def generate(self,registers):

        Rphys = QuantumRegister(self.n_phys,name='phys')
        Rtime = QuantumRegister(self.n_time,name='time')
        Rspecies = QuantumRegister(self.n_species,name='species')
        Rdir = QuantumRegister(self.n_dir,name='dir')
        Rsub = QuantumRegister(self.n_subnode,name='sub')
        Ra = QuantumRegister(self.n_ancila,name='a')
        
        QC = QuantumCircuit(Rphys,Rtime,Rspecies,Rdir,Rsub,Ra)
        
        registers.set_phys(Rphys[::-1])
        registers.set_time(Rtime[::-1])
        registers.set_species(Rspecies[::-1])
        registers.set_dir(Rdir[::-1])
        registers.set_sub(Rsub[::-1])
        registers.set_a(Ra[::-1])
        return QC

class Calculation_circuit_sampling():
    def __init__(self,Ndata):
        self.n_phys = Ndata[0]
        self.n_order = Ndata[1]
        self.n_species = Ndata[2]
        self.n_subnode = Ndata[3]
        self.n_ancila = Ndata[4]
        self.n_cla = Ndata[5]
    
    def generate(self,registers):

        Rphys = QuantumRegister(self.n_phys,name='phys')
        Rorder = QuantumRegister(self.n_order,name='order')
        Rspecies = QuantumRegister(self.n_species,name='species')
        Rsub = QuantumRegister(self.n_subnode,name='sub')
        Ra = QuantumRegister(self.n_ancila,name='a')
        cla = ClassicalRegister(self.n_cla,name='c')
        
        QC = QuantumCircuit(Rphys,Rorder,Rspecies,Rsub,Ra,cla)
        
        registers.set_phys(Rphys[::-1])
        registers.set_order(Rorder[::-1])
        registers.set_species(Rspecies[::-1])
        registers.set_sub(Rsub[::-1])
        registers.set_a(Ra[::-1])
        registers.set_cla(cla[::-1])
        return QC

class Registers():

    def set_phys(self,a):
        self.phys = [a[i] for i in range(len(a))][::-1]

    def set_order(self,a):
        self.order = [a[i] for i in range(len(a))][::-1]

    def set_species(self,a):
        self.species = [a[i] for i in range(len(a))][::-1]

    def set_sub(self,a):
        # self.sub = [a[i] for i in range(len(a))][::-1]
        self.sub = [a[i] for i in range(len(a))][::-1]
        
    def set_a(self,a):
        self.a = [a[i] for i in range(len(a))][::-1]

class Registers_m():

    def set_phys(self,a):
        self.phys = [a[i] for i in range(len(a))][::-1]

    def set_time(self,a):
        self.time = [a[i] for i in range(len(a))][::-1]
    
    def set_species(self,a):
        self.species = [a[i] for i in range(len(a))][::-1]

    def set_dir(self,a):
        self.dir = [a[i] for i in range(len(a))][::-1]

    def set_sub(self,a):
        # self.sub = [a[i] for i in range(len(a))][::-1]
        self.sub = [a[i] for i in range(len(a))][::-1]
        
    def set_a(self,a):
        self.a = [a[i] for i in range(len(a))][::-1]

    
class Registers_sampling():

    def set_phys(self,a):
        self.phys = [a[i] for i in range(len(a))][::-1]

    def set_order(self,a):
        self.order = [a[i] for i in range(len(a))][::-1]

    def set_species(self,a):
        self.species = [a[i] for i in range(len(a))][::-1]

    def set_sub(self,a):
        # self.sub = [a[i] for i in range(len(a))][::-1]
        self.sub = [a[i] for i in range(len(a))][::-1]
        
    def set_a(self,a):
        self.a = [a[i] for i in range(len(a))][::-1]

    def set_cla(self,a):
        self.cla = [a[i] for i in range(len(a))][::-1]