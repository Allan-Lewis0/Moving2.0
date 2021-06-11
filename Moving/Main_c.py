import numpy as np
import random as rd
import math
from Moving import Moving
import matplotlib.pyplot as plt

#实验可变参数
no_cell=50 #网格数
no_molecule_each_cell=30 #每网格分子数
no_molecule=1500#分子总数
jnis=50 #开始取样循环数
nloop=500 #总循环数
p=np.array([[0.0]*no_molecule for i in range(6)]) # p[0],p[1],p[2]为分子速度在x,y,z方向上的分量,p[3]为分子在y方向上的位置
vmean=np.array([[0.0]*no_molecule for i in range(3)]) #vmean[i,n],(0)(1)(2)为分子的信息速度在x,y,z方向上的分量
lcr=[0]*no_molecule #lcr[n]中存放新位置编码为n的分子的原始号码
yc=[0.0]*no_cell #yc[n]中存放第n网格中心的坐标
ic=np.array([[0]*no_cell for i in range(2)]) #ic[0,i]中存放第i网格中的分子数，ic[1,i]中存放第i网格中第1分子的序号-1
coll_remain=[0.0]*no_cell #coll_remain[i]中是第i网格中dtm时间间隔内剩余不足1的碰撞数
wall=[0.0]*7 #wall[i]是在下S平板累计的：（0）入射分子计数；（1）（4）入射和反射切向动量；（2）（5）入射和反射法向动量；（3）（6）入射和反射动能
field=np.array([[0.0]*no_cell for i in range(7)]) #流场中各量求和：（0,n）是第n网格中的分子数量累计；（1,n）（2,n）是第n网格中速度（*信息速度）1，2分量累计；（3,n）（4,n）（5,n）是第n网格中速度分量平方累计
op=[0.0]*8 #输出数据单元
#x,y,z坐标点集
x_values=[]
y_values=[]
z_values=[]
#缓冲区，用于选取点集中的部分点
x_values_std=[]
y_values_std=[]
z_values_std=[]


#Constant模块部分数据
pi = 3.1415926 #圆周率
rkn = 10.0 #KNUDSEN数的倒数
boltz = 1.3805e-23 #Boltzmann常数，k，（J/K）
viscosity = 2.117e-5 #氩的粘性系数，（Ns/m**2）
rmass = 6.63e-26 #氩的原子量，（kg）
t_ini = 273.0
pressure = 101325
fnd = pressure/(boltz*t_ini)
#以上三量为初始稳度，初始压力，初始数密度
u_wall = 100.0
# * u_wall = 1.0
#壁面速度
vm_ini = math.sqrt(2.0*boltz*t_ini/rmass) #初始相对速度
amda = viscosity*16.0/(5.0*math.sqrt(pi)*rmass*fnd*vm_ini) #硬球模型的分子平均自由程LAMBDA
y_length =4e-6  #rkn*amda
dtm = 0.23*amda/vm_ini
cell_height = (y_length)/float(no_cell) #网格高度
for n in range(0,no_cell):
    yc[n] = (n+0.5)*cell_height+(-2e-6) #第N网格中心的坐标
    coll_remain[n] = rd.random() #每网格的剩余碰撞数，初始为随机分数
print("初始化......")
print("壁面范围：0<=Y<={}".format(y_length))
print("每网格高度：",cell_height)
print("时间步长：",dtm)
print("第N网格中心的坐标：yc[]=",yc)
print("每网格的剩余碰撞数：coll_remain[]=",coll_remain)

#开始模拟分子运动过程
moving=Moving() #创建分子运动模型对象
iloop=0  #当前循环数
print("模拟分子运中......")
#置分子初速度与位置
moving.subc2(no_molecule,no_cell,no_molecule_each_cell,vm_ini,cell_height,yc,p,u_wall)
# for m in range(0,no_molecule):
#     print("p[3,{}]".format(m),p[3,m])
while(iloop<nloop): #当前循环数小于总循环数时，继续循环
    iloop=iloop+1
    if iloop%100==0:
        print("已运动{}个时间步长".format(iloop))
    for m in range(0,no_molecule):
        #记录X,Y,Z的坐标信息
        next_x=p[4,m]
        next_y=p[3,m]
        next_z=p[5,m]
        x_values.append(next_x)
        y_values.append(next_y)
        z_values.append(next_z)

    next_x1=x_values[(iloop-1)*no_molecule+1236]
    next_y1=y_values[(iloop-1)*no_molecule+1236]
    next_z1=z_values[(iloop-1)*no_molecule+1236]
    # print("iloop={},x=".format(iloop))
    # if iloop>1:
    #     print("dis=",math.sqrt((next_x-x_values_std[-1])**2+(next_y-y_values_std[-1])**2+(next_z-z_values_std[-1])**2))
    x_values_std.append(next_x1)
    y_values_std.append(next_y1)
    z_values_std.append(next_z1)
    # if iloop==1:
    #     for m in range(0,no_molecule):
    #         #记录X,Y,Z的坐标信息
    #         next_x=x_values[m]
    #         next_y=y_values[m]
    #         next_z=z_values[m]
    #         x_values_std.append(next_x)
    #         y_values_std.append(next_y)
    #         z_values_std.append(next_z)
    # if iloop==1000:
    #     for m in range(0,no_molecule):
    #         #记录X,Y,Z的坐标信息
    #         next_x=x_values[-m]  #x_values[-1]+p[0,m]*dtm
    #         next_y=y_values[-m]
    #         next_z=z_values[-m]  #z_values[-1]+p[2,m]*dtm
    #         x_values_std.append(next_x)
    #         y_values_std.append(next_y)
    #         z_values_std.append(next_z)


    #计算分子运动和在表面的反射
    # print("第{}次循环".format(iloop))
    moving.subc3(iloop,jnis,no_molecule,p,wall)
    #分子排序，编号
    moving.subc4(no_molecule,no_cell,cell_height,p,ic,lcr)
    #计算碰撞
    moving.subc5(no_cell,no_molecule,coll_remain,ic,lcr,p)
    if iloop>jnis: #当前循环数大于开始取样循环数时进行壁面和流场取样
        moving.subc6(no_cell,no_molecule,ic,lcr,p,field)

#当前循环数等于总循环数时进行取样和输出
moving.subc7(nloop,jnis,no_cell,yc,no_molecule_each_cell,wall,field,op)
print("{}个分子进行了{}个时间步长的运动，总共生成{}个点".format(no_molecule,iloop,len(x_values)))
print(len(x_values_std))

#绘制点的运动范围
ax = plt.subplot(111, projection='3d')  # 创建一个三维的绘图工程
# 将数据点分成三部分画，在颜色上有区分度
ax.scatter(x_values, y_values, z_values, s=0.5,c='y')  # 绘制数据点
ax.scatter(x_values[0], y_values[0], z_values[0], c='green')
ax.scatter(x_values[-1], y_values[-1], z_values[-1], c='red')
ax.set_zlabel('Z')  # 坐标轴
ax.set_ylabel('Y')
ax.set_xlabel('X')
# ax.view_init(0,0)
plt.show()