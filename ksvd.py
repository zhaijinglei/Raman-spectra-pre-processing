import numpy as np
from lmfit.models import (VoigtModel)  # pip install lmfit
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os
from scipy import interpolate
from make_dictionary import fit
from baselinewavelet import baselineWavelet
from scipy.sparse import csc_matrix, eye, diags
# from baselineCorrection import baseline_als
from scipy.stats import norm
def MMS(data):  # 最大值归一化
    return [float(i)  / (max(data)) for i in data]
# x=np.arange(0,1031,1)
# signal=np.zeros((100,1031))
# for i in range(100):
# 	g1=norm(loc = i*10, scale = 3.0)  
# 	signal[i]=g1.pdf(x)
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./test/p.csv', header=None)
######################################################################
# 硫高斯字典
# x=np.arange(0,1031,1)
# signal=np.zeros((5,1031))
# # g1=norm(loc =14, scale = 3.0)  
# g2 = norm(loc=84, scale=3.0)   
# g3 = norm(loc=150, scale=3.0)   
# g4 = norm(loc=404, scale=3.0)   
# g5 = norm(loc=177, scale=3.0)   
# g6 = norm(loc=366, scale=3.0)   
# signal[0]=g6.pdf(x)
# signal[1]=g2.pdf(x)
# signal[2]=g3.pdf(x)
# signal[3]=g4.pdf(x)
# signal[4]=g5.pdf(x)
# # signal[5]=g6.pdf(x)
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./test/p_s.csv', header=None)
######################################################################

######################################################################
# 对乙酰氨基酚高斯字典
# x=np.arange(0,1636,1)
# signal=np.zeros((20,1636))
# g1=norm(loc =214-65, scale = 3.0)
# g2 = norm(loc=328-65, scale=3.0)
# g3 = norm(loc=391-65, scale=3.0)
# g4 = norm(loc=465-65, scale=3.0)
# g5 = norm(loc=504-65, scale=3.0)
# g6 = norm(loc=651-65, scale=3.0)
# g7 = norm(loc=710-65, scale=3.0)
# g8 = norm(loc=797-65, scale=3.0)
# g9 = norm(loc=834-65, scale=3.0)
# g10 = norm(loc=857-65, scale=3.0)
# g11 = norm(loc=968-65, scale=3.0)
# g12 = norm(loc=1104-65, scale=3.0)
# g13 = norm(loc=1168-65, scale=3.0)
# g14 = norm(loc=1236-65, scale=3.0)
# g15 = norm(loc=1277-65, scale=3.0)
# g16 = norm(loc=1324-65, scale=3.0)
# g17 = norm(loc=1371-65, scale=3.0)
# g18 = norm(loc=1515-65, scale=3.0)
# g19= norm(loc=1561-65, scale=3.0)
# g20 = norm(loc=1648-65, scale=3.0)
#
# signal[0]=g1.pdf(x)
# signal[1]=g2.pdf(x)
# signal[2]=g3.pdf(x)
# signal[3]=g4.pdf(x)
# signal[4]=g5.pdf(x)
# signal[5]=g6.pdf(x)
# signal[6]=g7.pdf(x)
# signal[7]=g8.pdf(x)
# signal[8]=g9.pdf(x)
# signal[9]=g10.pdf(x)
# signal[10]=g11.pdf(x)
# signal[11]=g12.pdf(x)
# signal[12]=g13.pdf(x)
# signal[13]=g14.pdf(x)
# signal[14]=g15.pdf(x)
# signal[15]=g16.pdf(x)
# signal[16]=g17.pdf(x)
# signal[17]=g18.pdf(x)
# signal[18]=g19.pdf(x)
# signal[19]=g20.pdf(x)
#
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./test/p_dy.csv', header=None)
######################################################################

######################################################################
# 甲醇高斯字典
# x=np.loadtxt('./gauss/chunx.txt')
# signal=np.zeros((5,len(x)))
# g1=norm(loc =1024, scale = 3.0)
# g2 = norm(loc=1103, scale=3.0)
# g3 = norm(loc=1450, scale=3.0)
# g4 = norm(loc=1150, scale=3.0)
# g5 = norm(loc=197, scale=3.0)
#
#
# signal[0]=g1.pdf(x)
# signal[1]=g2.pdf(x)
# signal[2]=g3.pdf(x)
# signal[3]=g4.pdf(x)
# signal[4]=g5.pdf(x)
#
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./gauss/p_jiachun.csv', header=None)
######################################################################
# 乙醇高斯字典
# x=np.loadtxt('./gauss/chunx.txt')
# signal=np.zeros((7,len(x)))
# g1=norm(loc =1048, scale = 3.0)
# g2 = norm(loc=1094, scale=3.0)
# g3 = norm(loc=1450, scale=3.0)
# g4 = norm(loc=1271, scale=3.0)
# g5 = norm(loc=197, scale=3.0)
# g6 = norm(loc=427, scale=3.0)
# g7 = norm(loc=872, scale=3.0)
#
#
# signal[0]=g1.pdf(x)
# signal[1]=g2.pdf(x)
# signal[2]=g3.pdf(x)
# signal[3]=g4.pdf(x)
# signal[4]=g5.pdf(x)
# signal[5]=g6.pdf(x)
# signal[6]=g7.pdf(x)
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./gauss/p_yichun.csv', header=None)

######################################################################
# 乙醇高斯字典
# x=np.loadtxt('./gauss/chunx.txt')
# signal=np.zeros((13,len(x)))
# g1=norm(loc =970, scale = 3.0)
# g2 = norm(loc=1025, scale=3.0)
# g3 = norm(loc=1097, scale=3.0)
# g4 = norm(loc=1296, scale=3.0)
# g5 = norm(loc=1450, scale=3.0)
# g6 = norm(loc=1271, scale=3.0)
# g7 = norm(loc=205, scale=3.0)
# g8 = norm(loc=175, scale=3.0)
# g9 = norm(loc=318, scale=3.0)
# g10 = norm(loc=473, scale=3.0)
# g11 = norm(loc=762, scale=3.0)
# g12 = norm(loc=858, scale=3.0)
# g13 = norm(loc=885, scale=3.0)
#
#
# signal[0]=g1.pdf(x)
# signal[1]=g2.pdf(x)
# signal[2]=g3.pdf(x)
# signal[3]=g4.pdf(x)
# signal[4]=g5.pdf(x)
# signal[5]=g6.pdf(x)
# signal[6]=g7.pdf(x)
# signal[7]=g8.pdf(x)
# signal[8]=g9.pdf(x)
# signal[9]=g10.pdf(x)
# signal[10]=g11.pdf(x)
# signal[11]=g12.pdf(x)
# signal[12]=g13.pdf(x)
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./gauss/p_chun.csv', header=None)
# DDT高斯字典
# x=np.arange(200,2201,1)
# peak=np.loadtxt('.\DDT.txt')
# signal=np.zeros((len(peak),len(x)))
# for i in range(len(peak)):
# 	g = norm(loc=peak[i], scale=3.0)
# 	signal[i]=g.pdf(x)
#
#
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./gauss/DDT.csv', header=None)
# #萘高斯字典
# #x=np.loadtxt('./gauss/chunx.txt')
# x=np.arange(65,1701,1)
# signal=np.zeros((7,1636))
# g1=norm(loc =763, scale = 3.0)
# g2 = norm(loc=1381, scale=3.0)
# g3 = norm(loc=513, scale=3.0)
# g4 = norm(loc=1019, scale=3.0)
# g5 = norm(loc=1463, scale=3.0)
# g6 = norm(loc=1576, scale=3.0)
# g7 = norm(loc=1146, scale=3.0)
# signal[0]=g1.pdf(x)
# signal[1]=g2.pdf(x)
# signal[2]=g3.pdf(x)
# signal[3]=g4.pdf(x)
# signal[4]=g5.pdf(x)
# signal[5]=g6.pdf(x)
# signal[6]=g7.pdf(x)
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./gauss/p_nai.csv', header=None)

#环己烷高斯字典
#x=np.loadtxt('./gauss/chunx.txt')
# x=np.arange(65,1701,1)
# signal=np.zeros((7,1636))
# g1=norm(loc =801, scale = 3.0)
# g2 = norm(loc=1028, scale=3.0)
# g3 = norm(loc=1266, scale=3.0)
# g4 = norm(loc=1444, scale=3.0)
# g5 = norm(loc=1157, scale=3.0)
# g6 = norm(loc=426, scale=3.0)
# g7 = norm(loc=384, scale=3.0)
# signal[0]=g1.pdf(x)
# signal[1]=g2.pdf(x)
# signal[2]=g3.pdf(x)
# signal[3]=g4.pdf(x)
# signal[4]=g5.pdf(x)
# signal[5]=g6.pdf(x)
# signal[6]=g7.pdf(x)
# matrix2 = pd.DataFrame(signal)
# matrix2.to_csv('./gauss/p_c6.csv', header=None)
def pre_data(x, y):

	f1 = interpolate.interp1d(x, y, kind='linear')  # 插值模块
	x_pred = np.linspace(400, 1700, num=1301)
	y_pred = f1(x_pred)

	return x_pred, y_pred
y=np.loadtxt('./环己烷/c6-ann.txt')
# for j in range(len(x)):
#     x[i] = float(x[i])
path = "F:/laman/小拉曼/环己烷/"  # 待读取的文件夹
path_list = os.listdir(path)
path_list.sort()  # 对读取的路径进行排序
for i in range(len(path_list)):
    with open(path + path_list[i], 'r') as csvfile:
        reader = csv.reader(csvfile)
        column = [row[0] for row in reader]
        x = column
    for j in range(len(x)):
        x[j]=float(x[j])
x, y=pre_data(x, y)
np.savetxt('./环己烷/c6ann.txt',y)