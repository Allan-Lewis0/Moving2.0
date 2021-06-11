import random as rd
import math
import Constant
import csv
from Moving_Set import Wall

class Moving:

    #置分子初速度与位置
    def subc2(self,no_molecule,no_cell,no_molecule_each_cell,vm_ini,cell_height,yc,p,u_wall):
        '''置分子初速度与位置'''
        pi=Constant.pi
        for n in range(0,no_cell):
            for m in range(0,no_molecule_each_cell):
                iall=n*no_molecule_each_cell+m
                for i in range(0,3):
                    abc1=math.sqrt(-math.log(rd.random()))
                    abc2=2.0*pi*rd.random()
                    p[i,iall]=abc1*math.sin(abc2)*vm_ini
                #这一循环产生模拟分子的初始分子速度，即其初始宏观速度分量加上平衡气体中热运动速度分量，初始速度为0
                rand1=rd.random()
                if rand1>0.5:
                    p[3,iall]=yc[n]+0.5*rand1*cell_height
                else:
                    p[3,iall]=yc[n]-0.5*rand1*cell_height
                p[4,iall]=0.0
                p[5,iall]=0.0
                #给出模拟分子的在网格中的初始随机位置，这里考虑了减少方差原则

    #计算分子运动和在上下表面的反射和取样
    def subc3(self,iloop,jnis,no_molecule,p,wall):
        '''计算分子运动和在上下表面的反射和取样'''
        pi=Constant.pi
        dtm=Constant.dtm
        u_wall=Constant.u_wall
        vm_wall=Constant.vm_wall
        m_set = Wall()
        for m in range(0,no_molecule):
            dir=-p[1,m]/math.fabs(p[1,m])
            x1=p[4,m]
            y1=p[3,m]
            z1=p[5,m]
            x=p[4,m]+p[0,m]*dtm
            y=p[3,m]+p[1,m]*dtm  #y是模拟分子以p（1，m）运动经dtm时间后到达的位置
            z=p[5,m]+p[2,m]*dtm
            flag=m_set.isWall(x,y,z)
            if flag[0]==1:  #分子在壁面反射
                distance=m_set.Wall_distance(x1,y1,z1)
                if z1==0:
                    dtr=dtm-distance/math.fabs(p[1,m])
                    if iloop>=jnis:
                        p[0,m]=p[0,m]+u_wall/2.0
                        wall[0]=wall[0]+1.0
                        wall[1]=wall[1]+p[0,m]
                        wall[2]=wall[2]-p[1,m]
                        wall[3]=wall[3]+0.5*(p[0,m]**2+p[1,m]**2+p[2,m]**2)
                        #当当前循环数iloop大于开始取样循环数jnis时，进行了上平板取样（wall）,累计了入射分子贡献：（1）计数；（2）切向动量；（3）法向动量；（4）动能。
                    abc1=math.sqrt(-math.log(rd.random()))*vm_wall
                    abc2=2.0*pi*rd.random()
                    p[0,m]=abc1*math.sin(abc2)
                    p[2,m]=abc1*math.cos(abc2)
                    p[1,m]=-dir*math.sqrt(-math.log(rd.random()))*vm_wall
                    #计算了在上平板漫反射后的分子速度，反射后的IP速度在与平板相联系的坐标中为零
                    if iloop>=jnis:
                        wall[4]=wall[4]+p[0,m]
                        #在IP方法中wall(4)中累计的是反射分子的IP速度
                        wall[5]=wall[5]+p[1,m]
                        wall[6]=wall[6]+0.5*(p[0,m]**2+p[1,m]**2+p[2,m]**2)
                    #以上累计了反射分子贡献：（4）切向动量(为0)；（5）法向动量；（6）动能。
                    p[0,m]=p[0,m]-u_wall/2.0
                    if dir==-1:
                        p[3,m]=-flag[1]+p[1,m]*dtr
                        p[4,m]=x1+p[0,m]*dtr
                        p[5,m]=z1+p[2,m]*dtr
                    else:
                        p[3,m]=flag[1]+p[1,m]*dtr
                        p[4,m]=x1+p[0,m]*dtr
                        p[5,m]=p[2,m]*dtr
                    #以上从上平板反射的切向分子速度加上了上平板速度，切向和法向信息速度赋值。分子运到新位置。
                    continue
                else:
                    dtr=dtm-distance/math.sqrt(p[1,m]**2+p[2,m]**2)  #反射后剩余的分子运动时间    *
                    if iloop>=jnis:
                        p[0,m]=p[0,m]+u_wall/2.0
                        # **  vmean[0,m] = vmean[0,m]+u_wall/2.0
                        #这里将参考系转换为以上平板（以-u_wall/2.0运动）相联系
                        wall[0]=wall[0]+1.0
                        wall[1]=wall[1]+p[0,m]
                        # * wall[1] = wall[1]+vmean[0,m]
                        #在IP方法中u_wall[1]中累计的是入射分子的IP速度
                        wall[2]=wall[2]-p[1,m]
                        wall[3]=wall[3]+0.5*(p[0,m]**2+p[1,m]**2+p[2,m]**2)
                        #当当前循环数iloop大于开始取样循环数jnis时，进行了上平板取样（wall）,累计了入射分子贡献：（1）计数；（2）切向动量；（3）法向动量；（4）动能。
                    abc1=math.sqrt(-math.log(rd.random()))*vm_wall
                    abc2=2.0*pi*rd.random()
                    abc3=2.0*pi*rd.random()
                    abc4=2.0*pi*rd.random()
                    p[0,m]=abc1*math.sin(abc2)
                    p[2,m]=abc1*math.sin(abc3)
                    p[1,m]=abc1*math.cos(abc4)
                    #计算了在上平板漫反射后的分子速度，反射后的IP速度在与平板相联系的坐标中为零
                    if iloop>=jnis:
                        wall[4]=wall[4]+p[0,m]
                        #在IP方法中wall(4)中累计的是反射分子的IP速度
                        wall[5]=wall[5]+p[1,m]
                        wall[6]=wall[6]+0.5*(p[0,m]**2+p[1,m]**2+p[2,m]**2)
                    #以上累计了反射分子贡献：（4）切向动量(为0)；（5）法向动量；（6）动能。
                    k=y1/z1
                    p[0,m]=p[0,m]-u_wall/2.0
                    p[3,m]=k*flag[1]/math.sqrt(k**2+1)+p[1,m]*dtr
                    p[4,m]=x1+p[0,m]*dtr
                    p[5,m]=flag[1]/math.sqrt(k**2+1)+p[2,m]*dtr
                    continue
            # if flag[0]==2:
            #     dtr=(p[4,m]-(2e-6))/p[0,m] #反射后剩余的分子运动时间
            #     if iloop>=jnis:
            #         p[0,m]=p[0,m]+u_wall/2.0
            #         # **  vmean[0,m] = vmean[0,m]+u_wall/2.0
            #         #这里将参考系转换为以上平板（以-u_wall/2.0运动）相联系
            #         wall[0]=wall[0]+1.0
            #         wall[1]=wall[1]+p[0,m]
            #         # * wall[1] = wall[1]+vmean[0,m]
            #         #在IP方法中u_wall[1]中累计的是入射分子的IP速度
            #         wall[2]=wall[2]-p[1,m]
            #         wall[3]=wall[3]+0.5*(p[0,m]**2+p[1,m]**2+p[2,m]**2)
            #         #当当前循环数iloop大于开始取样循环数jnis时，进行了上平板取样（wall）,累计了入射分子贡献：（1）计数；（2）切向动量；（3）法向动量；（4）动能。
            #     abc1=math.sqrt(-math.log(rd.random()))*vm_wall
            #     abc2=2.0*pi*rd.random()
            #     p[1,m]=abc1*math.sin(abc2)
            #     p[2,m]=abc1*math.cos(abc2)
            #     p[0,m]=math.sqrt(-math.log(rd.random()))*vm_wall
            #     #计算了在上平板漫反射后的分子速度，反射后的IP速度在与平板相联系的坐标中为零
            #     if iloop>=jnis:
            #         wall[4]=wall[4]+p[0,m]
            #         # * wall[4] = wall[4] + vmean[0, m]
            #         #在IP方法中wall(4)中累计的是反射分子的IP速度
            #         wall[5]=wall[5]+p[1,m]
            #         wall[6]=wall[6]+0.5*(p[0,m]**2+p[1,m]**2+p[2,m]**2)
            #     #以上累计了反射分子贡献：（4）切向动量(为0)；（5）法向动量；（6）动能。
            #     p[0,m]=p[0,m]-u_wall/2.0
            #     # **  vmean[0,m] = -u_wall/2.0
            #     # ** vmean[1, m] = 0.0
            #     p[3,m]=y1+p[1,m]*dtr
            #     p[4,m]=(2e-6)+p[0,m]*dtr
            #     p[5,m]=z1+p[2,m]*dtr
            #     #以上从上平板反射的切向分子速度加上了上平板速度，切向和法向信息速度赋值。分子运到新位置。
            #     continue
            else:
                p[3,m]=y
                p[4,m]=x
                p[5,m]=z
            continue

    #原始号码为M的分子根据其运动后所在网格重新编号为K.M存于LCR(K)中。
    def subc4(self,no_molecule,no_cell,cell_height,p,ic,lcr):
        '''原始号码为M的分子根据其运动后所在网格重新编号为K.M存于LCR(K)中。'''

        for icell in range(0,no_cell):
            ic[0,icell]=0
        for m in range(0,no_molecule):
            ncell=int((p[3,m]+2e-6)/cell_height)
            if ncell>=no_cell:
                ncell=no_cell-1  #若写成ncell = no_cell，则ic[0,50]超出界限，最大只到ic[0,49]
                #NCELL是M分子所在的新的网格的编号。
            ic[0,ncell]=ic[0,ncell]+1
        #这里数好了所有NO_CELL个网格中的分子数IC[0,I]，为得到IC[1,I]做准备。
        ic2=0
        for icell in range(0,no_cell):
            ic[1,icell]=ic2
            #现在IC[1,I]已经是第I网格第1个分子的序号
            ic2=ic2+ic[0,icell]
            ic[0,icell]=0
        for m in range(0,no_molecule):
            ncell=int((p[3,m]+2e-6)/cell_height)
            if ncell>=no_cell:
                ncell=no_cell-1
            ic[0,ncell]=ic[0,ncell]+1
            k=int(ic[1,ncell]-1+ic[0,ncell])
            lcr[k]=m
            #此循环又数出了IC[0,I]，并得到了原始号码为M的分子根据其新位置的编号K。将M置于LCR[K]中。

    #计算碰撞
    def subc5(self,no_cell,no_molecule,coll_remain,ic,lcr,p):
        '''计算碰撞'''
        vrc=[0.0]*3
        vccm=[0.0]*3
        pi=Constant.pi
        eta51=Constant.eta51
        area=Constant.area
        cell_height=Constant.cell_height
        tkom=Constant.tkom
        dtm=Constant.dtm
        vrmax_eta51=Constant.vrmax_eta51

        for m in range(0,no_cell):
            if ic[0,m]<2:
                continue
                #如果M网格中无分子或仅有一个分子，跳过碰撞
            vaver=0.0  #vaver是网格中CR**ETA51累计
            icontr=ic[0,m]
            if ic[0,m]==2:
                icontr=1
            if ic[0,m]==3:
                icontr=2
            #icontr是需要随机选中的分子对数 ，以计算网格中DTM内的碰撞数
            for icoll in range(0,icontr):
                k=int(rd.random()*ic[0,m]+ic[1,m])
                l1=lcr[k]  #在M网格在随机选中一个分子l1
                j=int(rd.random()*ic[0,m]+ic[1,m])
                while (j==k):
                    j=int(rd.random()*ic[0,m]+ic[1,m])
                l2=lcr[j]  #在M网格在随机选中另一个分子l2

                for i in range(0,3):
                    vrc[i]=p[i,l1]-p[i,l2]
                vr=math.sqrt(vrc[0]**2+vrc[1]**2+vrc[2]**2)
                vaver=vaver+vr**eta51
                #vr是L1和L2分子的相对速度，所得的VR**eta51乘以tkom即得s(cr)
            vaver=vaver/icontr
            cnoic=ic[0,m]*(ic[0,m]/(area*cell_height))*tkom*vaver*dtm/2.0+coll_remain[m]
            #大括弧为流场数密度n，tkom*vaver为S[CR]的平均值，cnoic即dtm时间内每网格内发生的碰撞数Nnsf
            ncoll=int(cnoic)
            #ncoll是考虑到上一dtm剩余碰撞数后再取整数为实际应发生的碰撞数
            # print("ic[0,{}]=".format(m),ic[0,m])
            # print("ncoll={}".format(ncoll))
            coll_remain[m]=cnoic-ncoll
            if ncoll<1:
                continue
            nacoll=0  #nacoll用来统计实际发生的碰撞数
            while nacoll<ncoll:
                k=int(rd.random()*ic[0,m]+ic[1,m])
                l1=lcr[k]
                j=int(rd.random()*ic[0,m]+ic[1,m])
                while j==k:
                    j=int(rd.random()*ic[0,m]+ic[1,m])
                l2=lcr[j]
                # 从网格中随机选出两个分子l1和l2
                for i in range(0,3):
                    vrc[i]=p[i,l1]-p[i,l2]
                vr=math.sqrt(vrc[0]**2+vrc[1]**2+vrc[2]**2)
                # vr为碰撞对的相对速度
                if vr**eta51>vrmax_eta51:
                    vrmax_eta51=vr**eta51
                a=vr**eta51/vrmax_eta51
                b=rd.random()
                if a<b:
                    continue
                #利用取舍法判断分子对是否中选，vr**eta51在eta51合适的赋值时可以包括VHS,VSS模型
                c=1.0-2.0*rd.random()
                d=math.sqrt(1.0-c*c)
                vrc[0]=c*vr
                t=2.0*pi*rd.random()
                vrc[1]=d*vr*math.cos(t)
                vrc[2]=d*vr*math.sin(t)
                # vrc为碰撞后的相对速度的三分量，适用于HS，VHS模型
                for k in range(0,3):
                    vccm[k]=0.5*(p[k,l1]+p[k,l2])
                    p[k,l1]=vccm[k]+vrc[k]*0.5
                    p[k,l2]=vccm[k]+vrc[k]*0.5
                nacoll=nacoll+1

    #流场参数求和
    def subc6(self,no_cell,no_molecule,ic,lcr,p,field):
        '''流场参数求和'''

        for icell in range(0,no_cell):
            for molecule in range(0,ic[0,icell]):
                k=ic[1,icell]+molecule
                l=lcr[k]
                field[0,icell]=field[0,icell]+1
                field[1,icell]=field[1,icell]+p[0,l]
                # * field[1,icell] = field[1,icell]+vmean[0,l]
                # IP方法中累计的是IP速度
                field[2,icell]=field[2,icell]+p[1,l]
                # * field[2,icell] = field[2,icell]+vmean[1,l]
                for i in range(0,3):
                    field[3+i,icell]=field[3+i,icell]+p[i,l]**2

    #输出相关参数
    def subc7(self,nloop,jnis,no_cell,yc,no_molecule_each_cell,wall,field,op):
        '''输出表面参量，记入surf_coue.dat;输出流场参量，记入u_coue.dat'''
        rkn=Constant.rkn
        u_wall=Constant.u_wall
        vm_ini=Constant.vm_ini
        dtm=Constant.dtm
        fnd=Constant.fnd
        area=Constant.area
        pi=Constant.pi
        rmass=Constant.rmass
        boltz=Constant.boltz
        t_ini=Constant.t_ini

        #以下输出壁面数据
        skn=1.0/rkn
        sma=math.sqrt(2.0/1.6670)*u_wall/vm_ini
        op[0]=wall[0]  #[0]分子计数，样本大小，（nloop-jnis）*dtm=t时间内打到上平板的分子数
        a=(nloop-jnis)*dtm*(fnd*area*vm_ini)
        op[1]=wall[0]/a  #[2]用n*vm无量纲化的数目通量wall(1)/t*A
        op[2]=wall[2]/(vm_ini*a)
        op[3]=wall[5]/(vm_ini*a)  #[2],[3]用n*vm**2无量纲化的入射和反射法向动量通量
        op[4]=2.0*math.sqrt(pi)*wall[1]/(u_wall*a)
        op[5]=2.0*math.sqrt(pi)*wall[4]/(u_wall*a)  #[4],[5]用TAUfm= m*n*U*vm/2*pi**0.5无量纲化的入射和反射法向动量通量。注意在IP方法中wall[1]与wall[4]已是入射和反射IP速度的累计
        op[6]=wall[3]/(a*vm_ini**2)
        op[7]=wall[6]/(a*vm_ini**2)  #[6],[7]用n*vm**3无量纲化的入射和反射法向能量通量
        drag=rmass*wall[1]/((nloop-jnis)*dtm*area)
        f_wall=open("E:\\pythonProject1\\Moving\\wall_data.csv","a",newline='')
        csv_writer1=csv.writer(f_wall,dialect="excel")
        csv_writer1.writerow(op)
        f_wall.close()

        #以下输出流场数据
        for n in range(0,no_cell):
            op[0]=field[0,n]  #网格中分子计数，样本大小
            op[1]=op[0]/((nloop-jnis)*no_molecule_each_cell)  #无纲量密度
            if op[0]==0:
                for i in range(2,8):
                    op[i]=0
            else:
                op[2]=field[1,n]/op[0]
                op[3]=field[2,n]/op[0]  #[2],[3]为平均切向，法向速度，* 平均切向，法向信息速度
                op[4]=rmass*(field[3,n]/op[0]-op[2]**2)/boltz/t_ini
                op[5]=rmass*(field[4,n]/op[0]-op[3]**2)/boltz/t_ini
                op[6]=rmass*field[5,n]/op[0]/boltz/t_ini
                #[4],[5],[6] X,Y,Z方向温度
                op[7]=(op[4]+op[5]+op[6])/3
            f_field=open("E:\\pythonProject1\\Moving\\field.csv","a",newline='')
            csv_writer3=csv.writer(f_field,dialect="excel")
            csv_writer3.writerow(op)
            f_field.close()
