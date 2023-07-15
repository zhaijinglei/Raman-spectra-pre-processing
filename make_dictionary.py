import numpy as np
from lmfit.models import (VoigtModel)  # pip install lmfit
from lmfit.models import (GaussianModel)
import matplotlib.pyplot as plt
import pandas as pd
from baselinewavelet import baselineWavelet


def iteration(data):  # 迭代补峰算法
    p = len(data)
    data[0] = 0
    #data[1999] = 0
    for i in range(20):
        data[12+i]=0
    # for i in range(400):
    #     data[i] = 0
    data[-1] = 0
    a = 1#=判断值，初始为1，检测到峰值则为0
    peakpoint = 1#初始化峰值点、左右边界点
    boundl = 0
    boundr = 0
    snn=0.001*max(data)#大拉曼强度值高，有可能陷入无限循环，截止条件需要调整
    mm = max(data)#找残差的最高点，在最高点加入一个峰
    for i in range(p - 1):
        if data[i] == mm and mm>snn:#默认峰值大于0.01倍的最大值，定义残差小于0.01倍的最大值时，退出检峰循环
            peakpoint = i#记录峰值点位置
            a = 0#检测到峰值点
            print(peakpoint )
            break
    Di = np.diff(data)#生成原始信号的导数矩阵
    for j in range(peakpoint):
        if (peakpoint + j) < p - 2:#右边界不超过光谱范围
            if (Di[peakpoint + j] < Di[peakpoint + j + 1]) and Di[peakpoint + j + 1] > -2e-190:#斜率绝对值从大变小，且绝对值小于一定范围，判断为右边界
                boundr = j + 1
                break
    for j in range(peakpoint):
        if (Di[peakpoint - j] > Di[peakpoint - j - 1]) and Di[peakpoint - j - 1] < 2e-190:#左边界思路亦然
            boundl = j + 1
            break
    for i in range(peakpoint - boundl):#将非峰区域置为0
        data[i] = 0
    for i in range(peakpoint):
        if Di[i] < 0:
            data[i + 1] = 0
    for i in range(p - peakpoint - boundr):#将非峰区域置为0
        data[peakpoint + boundr + i] = 0
    return data, boundr, boundl, a, peakpoint

def Voigt(x,y):
    x_shifted = x - x.min()  # Shifting to 0
    y_shifted = y - y.min()  # Shifting to 0
    mod = VoigtModel()  # Setting model type
    pars = mod.guess(y_shifted, x=x_shifted)  # Estimating fit
    out = mod.fit(y_shifted, pars, x=x_shifted)  # Fitting fit
    fwhm = out.params['fwhm'].value
    return out.best_fit,fwhm

def MMS(data):  # 最大值归一化
    return [float(i)  / (max(data)) for i in data]

def fit(x,data1):
    n = 10000  # 定义最大计算次数
    nn = len(data1)
    data1[1]=0
    peak = x - x
    save = peak
    # plt.figure()
    # plt.plot(x,data1, '-k')
   
    for i in range(nn):#去除过小的杂乱信号
        if data1[i] <20:
            data1[i] = 0
            
    # 判断光谱中峰的个数
    pea = peak
    for i in range(n):

        data = data1 - pea
        peak2, right, left, stop, peakp = iteration(data)
        if stop == 1:
            num = i
            break

        pea = pea + peak2

    p = np.arange(0, num, 1)
    matrix = np.ones((num, nn))
    # 拟合洛伦兹峰，并保存其和峰值点位
    for i in range(num):
        data = data1 - peak
        peak2, right, left, stop, peakp = iteration(data)
        x, y = np.array(x), np.array(peak2)
        y1 = y - y

        u, v = Voigt(x,y)  # 洛伦兹峰拟合
        # if v > 100:
        #     break
        y = y - u
        for j in range(len(x)):
            if y[j] <= 0:
                y[j] = 0
        y1 = y1 + u
        save = save + y1
        # 保存归一化标准峰
        matrix[i] = y1
        #plt.plot(x, y1, '-b')
        peak = peak + peak2
        p[i] = peakp
    return save, matrix,p


# if __name__ == '__main__':
#
#     y2 = np.loadtxt('test2.txt')
#     y2 = baselineWavelet(y2)
#     y2[y2 < 0] = 0
#     y2 = np.array(y2)
#     y22,matrix2 = fit(y2)
#
#
#     matrix2 = pd.DataFrame(matrix2)
#     matrix2.to_csv('./peak2.csv',header=None)
