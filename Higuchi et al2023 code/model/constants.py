import numpy as np

xmax=1*10**3#m
xmin=0

ymax=1*10**3#m
ymin=0

zmax=1*10**3#m
zmin=0

vxmax=1*10**8#m/s
vxmin=0

vymax=1*10**8#m/s
vymin=0

vzmax=1*10**8#m/s
vzmin=0

#各定数(SI単位系)MLTC単位系
e = 1.602*10**-19#*10**(-19)#電気素量
me = 9.10938*10**(-31)#kg質量
mi = 1000*me
c =299792458#m/s

#密度正規分布の定数
sigma =5*10**-3#分散
C1 = 1/(2*np.pi*sigma**2)
C2 = 2*sigma**2

#マクスウェル分布の定数

kT = 3*10**3*e#2*10**2*e#eV
C3=((me/(2*np.pi*kT))**(1/2))
C4=me/(2*kT)

#格子条件
Δt = 10**(-10)
Δx = c*Δt
Δy = Δx
Δz = Δx
Δvx = 10**8
Δvy = Δvx
Δvz = Δvx