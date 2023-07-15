#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'拟合过程中应当注意，坐标的顺序'

__author__ = 'alluring'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate

def fit(x,y,fit_num):#fit_num表示自由度
    params = np.polyfit(x, y, fit_num)  # 用5项多项式来拟合
    funcs = np.poly1d(params)  # funcs为拟合函数
    ypre = funcs(x)  # 用拟合函数和x值来预期y值
    y=y-ypre
    y[y<0]=0
    # np.savetxt('duoxiangshi.txt',y)
    plt.plot(y,'-r')
    plt.show()
    return y
def pre_data(x, y):
	f1 = interpolate.interp1d(x, y, kind='linear')  # 插值模块
	x_pred = np.linspace(400, 1700, num=1301)
	y_pred = f1(x_pred)

	return x_pred, y_pred

import csv
data_name = 'F:/laman/小拉曼/环己烷/c6-1s-1-50mw-10.csv'

with open(data_name, 'r') as csvfile:
    reader = csv.reader(csvfile)
    data = [row[1] for row in reader]  #
for i in range(len(data)):
    data[i] = float(data[i])
with open(data_name, 'r') as csvfile:
    reader = csv.reader(csvfile)
    xp = [row[0] for row in reader]  #
for i in range(len(xp)):
    xp[i] = float(xp[i])

x,y=pre_data(xp,data)
# x=np.arange(200,2201,1)
# y=np.loadtxt('./test/test.txt')
y=fit(x,y,7)
# df = pd.read_csv('F:/laman/小拉曼/对乙酰氨基酚/dy-1s-1-50mw-1.csv')
# data = df.values
# data = data.T
# x = data[0]
# y = data[1]
np.savetxt('./环己烷/环己烷 -PF.txt',y)


