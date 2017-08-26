#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
import copy
# from sympy import latex


# h为单元的厚度
# u为泊松比
h, u = 1, 0

# E为弹性模量
E = Symbol('E')

#节点荷载的列阵
Rx1, Rx2, Rx4,Ry4,Ry5,Ry6 = symbols('Rx1, Rx2, Rx4,Ry4,Ry5,Ry6')
Q = Matrix([Rx1,-10, Rx2, 0, 0, 0, Rx4,Ry4,0,Ry5,0,Ry6])


#节点位移的列阵
v1, v2, u3, v3, u5, u6 = symbols('v1 v2 u3 v3 u5 u6')

s = [0, v1 , 0, v2 , u3 , v3 , 0, 0, u5 , 0, u6 , 0]

def DimensionalityReduction(KK):
    # 需要删除矩阵的行列（找出列表中为零的序号）
    l = []
    nx = 0
    for i in s:
        nx += 1
        if i == 0:
            l.append(nx - 1)
    # print l
    # 矩阵降维
    n = 0
    for i in l:
        KK.col_del(i - n)
        KK.row_del(i - n)
        Q.row_del(i - n)
        n += 1

    KK = KK.col_insert(KK.shape[1], Q)
    d = KK.rref()[0].col(-1)

    # 有确定值的位移列阵
    for i in l:
        d = d.row_insert(i,Matrix([[0]]))
    # pprint(d)
    return d

class Plane2D(object):
    def __init__(self,xl,xm,xn,yl,ym,yn,l,m,n,KK):
        self.l = l
        self.m = m
        self.n = n

        self.xl = xl
        self.xm = xm
        self.xn = xn

        self.yl = yl
        self.ym = ym
        self.yn = yn

        self.KK = KK
        """坐标的矩阵"""
        self.A = Matrix([[1,xl,yl], [1, xm, ym], [1, xn,yn]])
        self.Area = 0.5 * self.A.det()

        """A的代数余子式"""
        Al = []
        for i in range(0, 3):
            for j in range(1, 3):
                Al.append(self.calCofactor(self.A, i, j))

        self.B = 1 / (2 * self.Area) * Matrix([[Al[0], 0, Al[2], 0, Al[4], 0], [0, Al[1], 0, Al[3], 0, Al[5]], \
                                               [Al[1], Al[0], Al[3], Al[2], Al[5], Al[4]]])

        self.D = (E/(1-u**2)) * Matrix([[1, u, 0], [u, 1, 0], [0, 0, (1-u)*0.5]])

    """计算代数余子式"""
    def calCofactor(self,A,i,j):
        b = copy.deepcopy(A)
        b.row_del(i)
        b.col_del(j)
        if (i+j)%2 == 0:
            return b.det()
        else:
            return -b.det()

    """计算总刚"""
    def krs(self):

        Area = self.Area
        A = self.A
        KK = self.KK
        B = self.B

        Bl = B[ : , 0 : 2]
        Bm = B[ : , 2 : 4]
        Bn = B[ : , 4 : 6]
        # pprint (B)
        BList = {self.l:Bl, self.m:Bm, self.n:Bn}

        for r in BList:
            for s in BList:
                krs = BList[r].T * self.D * BList[s] * h * Area
                for i in range(2):
                    for j in range(2):
                        KK[2 * r + i - 2, 2 * s + j - 2] += krs[i, j]
        return KK

	# """计算单元应力"""
    def Stress(self,d):

        l, m, n = self.l, self.m, self.n
        A = Matrix([[1, self.xl, self.yl], [1, self.xm, self.ym], [1, self.xn, self.yn]])
        B = self.B
        D = (E / (1 - u ** 2)) * Matrix([[1, u, 0], [u, 1, 0], [0, 0, (1 - u) * 0.5]])

        #k单元的位移
        k = Matrix([0])
        for i in [l,m,n]:
            k = k.row_insert(-1,d.row(2*i-2))
            k = k.row_insert(-1, d.row(2 * i - 1))
        k.row_del(-1)
        # pprint(k)
        pprint (D*B*k)

if __name__ == "__main__":
	# 单元的总刚
    KK = zeros(12, 12)

    # 荷载Q
    # xl,xm,xn,yl,ym,yn,l,m,n
    xy = [[0, 0, 1, 2, 1, 1],[0, 1, 1, 1, 0, 1],[0, 0, 1, 1, 0, 0],\
          [1, 1, 2, 1, 0, 0]]
		  
	# 单元的结点号
    jdh = [[1, 2, 3],[2, 5, 3],[2, 4, 5],[3, 5, 6]]

    for i in range(4):
        KK = Plane2D(xy[i][0],xy[i][1],xy[i][2],xy[i][3],xy[i][4],xy[i][5],\
                jdh[i][0],jdh[i][1],jdh[i][2],KK).krs()
    # pprint (KK*1/E)
    # pprint (KK[1])
	
    d = DimensionalityReduction(KK)
    # pprint (d*E)
    for i in range(4):
        Plane2D(xy[i][0],xy[i][1],xy[i][2],xy[i][3],xy[i][4],xy[i][5],\
                jdh[i][0],jdh[i][1],jdh[i][2],KK).Stress(d)

    # pprint(Q)
