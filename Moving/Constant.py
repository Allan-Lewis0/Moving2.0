import math

no_molecule=1500#分子总数
no_cell=50 #网格数
pi = 3.1415926 #圆周率
rkn = 10.0 #KNUDSEN数的倒数
boltz = 1.3805e-23 #Boltzmann常数，k，（J/K）
viscosity = 2.117e-5 #氩的粘性系数，（Ns/m**2）
tref = 273.0 #参考温度度
rmass = 6.63e-26 #氩的原子量，（kg）
diaref = 3.659e-10 #氩的HS分子直径，（m），如用VHS模型，取diaref=4.17e-10
power = 999999.9 #即ETA，逆幂率中的幂次，这里之值相应硬球模型，如用VHS模型，power=7.452
eta51 = (power-5.0)/(power-1.0) #逆幂率分子的S(CR)，依赖于CR的幂次
ksai = 2.0/(power-1.0)
amr = rmass/2.0 #折合质量
cxsref = pi*diaref**2
vhs_core = 1.0
tkom = cxsref*(2.0*boltz*tref/amr)**ksai/vhs_core
t_ini = 273.0
pressure = 101325
fnd = pressure/(boltz*t_ini)
#以上三量为初始稳度，初始压力，初始数密度
t_wall = 273.0
u_wall = 100.0
# * u_wall = 1.0
#以上二量为壁面温度与速度
vm_ini = math.sqrt(2.0*boltz*t_ini/rmass)  #最概然速率，又称“最可几速率”，当气体处于热力学平衡态，分子符合麦克斯韦速率分布，与麦克斯韦速率分布f(v)的极大值对应。
vm_wall = math.sqrt(2.0*boltz*t_wall/rmass)
vrm_ini = 2.0*math.sqrt(2.0/pi)*vm_ini

#以上三量为初始的，壁面的最可几速度和初始的相对速度
vrmax_eta51 = 2.0*vm_ini**eta51 #S(CR)的初始最大值SMAX
amda = viscosity*16.0/(5.0*math.sqrt(pi)*rmass*fnd*vm_ini) #硬球模型的分子平均自由程LAMBDA：一个气体分子在连续两次碰撞之间可能通过的各段自由程的平均值

y_length = rkn*amda
dtm = 0.23*amda/vm_ini
#以上为流场范围和时间步长
area = no_molecule/(fnd*2e-6) #area是相对于模拟分子数的，一维问题中我们观察流场的代表截面积
cell_height = (2e-6)/float(no_cell) #网格高度

