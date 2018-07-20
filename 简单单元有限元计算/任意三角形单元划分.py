#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
import copy
import matplotlib.pyplot as plt
# from Tkinter import *
import argparse

h, u = 1, 0
E = Symbol('E')

#节点荷载的列阵
Rx1, Rx2, Rx4,Ry4,Ry5,Ry6 = symbols('Rx1, Rx2, Rx4,Ry4,Ry5,Ry6')

#节点位移的列阵
v1, v2, u3, v3, u5, u6 = symbols('v1 v2 u3 v3 u5 u6')

s = []
def JiangWei(KK, n, Q):
    # l需要删除矩阵的行列（找出列表中为零的序号）
    #nn为划分的数

    l = []

    for i in xrange(n + 1):
        node = TTy(i)
        if i == n:
            l.append(2 * (node - 1))
            for j in xrange(node, node + n + 1):
                l.append(2 * j - 1)
        else:
            l.append(2 * (node - 1))
    # print l
    # print len(l)

    # 矩阵降维
    p = 0
    for i in l:
        KK.col_del(i - p)
        KK.row_del(i - p)
        p += 1
    # pprint (KK/E)
    # print KK.shape
    for i in xrange(KK.shape[0]-1):
        Q = Q.row_insert(1, Matrix([[0]]))

    # 解方程、位移的解
    # d为方程的解
    # pprint (KK)
    KK = KK.col_insert(KK.shape[1], Q)
    d = KK.rref()[0].col(-1)

    # 有确定值的位移列阵(完整的）
    for i in l:
        d = d.row_insert(i,Matrix([[0]]))
    # print d
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
        for i in xrange(0, 3):
            for j in xrange(1, 3):
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

        BList = {self.l:Bl, self.m:Bm, self.n:Bn}

        for r in BList:
            for s in BList:
                krs = BList[r].T * self.D * BList[s] * h * Area
                for i in xrange(2):
                    for j in xrange(2):
                        KK[2 * r + i - 2, 2 * s + j - 2] += krs[i, j]
        return KK

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
        # t.insert(END,'sj')
        # return k

# def printStress(self, k):
#     pprint (self.D*self.B*k)

def TTy(n, b=1, c=0):
    if n == c:
        return b+n
    else:
        return TTy(n, b=b+c, c=c+1)

def TTx(n, b=1, c=0):
    if n == c:
        return b+n
    else:
        return TTx(n, b=b+c+1, c=c+1)


#def TTx(n):
#   if n == 0:
#        return 1
#    else:
#        return n+1+TTx(n-1)

def dividedUnit(n):
    """
    num是点的节点号码
    n是划分的单元数
    d是一小块的边长
    """

    d = 2.0/n
    y, x = 2.0, 0.0

    xylist = [0]

    #画点、点的注释
    for i in xrange(n+1):
        y = 2 - i*d
        x = 0.0
        for j in xrange(i+1):
            x = j * d
            xylist.append((x, y))

    # print len(xylist)
    # print xylist
    #为单元划分
    unitlist = []
    jdhlist = []
    for i in xrange(1, n + 1):
        xStart = TTy(i)
        xEnd = TTy(i)+i

        sStart = TTy(i - 1)
        flag = 1

        while not xStart == xEnd:

            if flag % 2 == 0:
                unitlist.append([xylist[sStart], xylist[xStart], xylist[sStart + 1]])
                jdhlist.append([sStart, xStart, sStart + 1])
                sStart += 1
                flag += 1

            else:
                unitlist.append([xylist[sStart], xylist[xStart], xylist[xStart + 1]])
                jdhlist.append([sStart, xStart, xStart + 1])
                xStart += 1
                flag += 1

    return unitlist, jdhlist

def drawingPlot(n, xy):
    num = 1
    d = 2.0 / n
    # print 'd=%.2f' %d
    y, x = 2.0, 0.0

    # 画点、点的注释
    for i in xrange(n + 1):
        y = 2 - i * d
        x = 0.0
        for j in xrange(i + 1):
            x = j * d
            plt.plot(x, y, 'ro')
            plt.text(x, y, num, color='blue', fontsize=10)
            num+=1

    for i in xrange(0,n):
        plt.plot([0,2-i*d], [2-i*d,0],'r')
        plt.plot([i*d,i*d],[0,2-(i)*d],'r')
        plt.plot([0,d*(i+1)],[2-(i+1)*d,2-(i+1)*d],'r')

    # print len(xy)
    num=1
    for i in xrange(len(xy)):
        x = (xy[i][0][0]+xy[i][1][0]+xy[i][2][0])/3.0
        y = (xy[i][0][1]+xy[i][1][1]+xy[i][2][1])/3.0
        plt.text(x,y,num,color='green', fontsize=15)
        num+=1

    plt.axis([0, 2, 0, 2])
    plt.grid(color='b', linewidth=0.5, linestyle='--')

    plt.show()

def printStress(xy, jdh,d):
    for i in xrange(len(jdh)):
        Plane2D(xy[i][0][0],xy[i][1][0],xy[i][2][0],xy[i][0][1],xy[i][1][1],xy[i][2][1],\
                jdh[i][0],jdh[i][1],jdh[i][2],KK).Stress(d)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", action="store", help="Divided unit", dest='n', type=int )

    parser.add_argument("-k",action="store_false", help="zong gang")
    parser.add_argument("-s", action="store_false", help="printing stress")
    parser.add_argument("-d", action="store_false", help="drawing picture")
    args = parser.parse_args()

    node = 2*(TTy(args.n)+args.n)
    KK = zeros(node, node)
    xy, jdh = dividedUnit(args.n)
    Q = Matrix([-10])
    # xl,xm,xn,yl,ym,yn,l,m,n

    for i in xrange(len(jdh)):
        KK = Plane2D(xy[i][0][0],xy[i][1][0],xy[i][2][0],xy[i][0][1],xy[i][1][1],xy[i][2][1],\
                     jdh[i][0],jdh[i][1],jdh[i][2],KK).krs()
    if not args.k:
        # print "------------总刚------------"
        pprint (KK)
    d = JiangWei(KK, args.n, Q)
    if not args.s:
        # print "------------应力------------"
        printStress(xy,jdh,d)
    if not args.d:
        # d = JiangWei(KK, args.n, Q)
        drawingPlot(args.n,xy)

    # pprint (KK*1/E)

    # pprint (d*E)
    # print '----------------------------------------'
    # drawingPlot(n, xy)

