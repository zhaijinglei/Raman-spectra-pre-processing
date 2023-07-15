import numpy as np
from lmfit.models import (VoigtModel)  # pip install lmfit

import pandas as pd

from whittaker_smooth import whittaker_smooth
import matplotlib.pyplot as plt
import matplotlib,os,csv
from scipy import interpolate
from make_dictionary import fit
from baselinewavelet import baselineWavelet
from scipy.sparse import csc_matrix, eye, diags
# from baselineCorrection import baseline_als
from scipy.stats import norm
from scipy.signal import savgol_filter

def pre_data(x, y):
	f1 = interpolate.interp1d(x, y, kind='linear')  # 插值模块
	x_pred = np.linspace(65, 1700, num=1636)
	y_pred = f1(x_pred)

	return x_pred, y_pred

def cosine_similarity(x, y):#余弦相似度
    num = x.dot(y.T)
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    n=num / denom
    if n>1:n=1.0
    return n


def judge(data):
	y=data
	path = "./lib-spectra/"  # 数据库中存储物质的光谱
	path_list = os.listdir(path)
	path_list.sort()  # 对读取的路径进行排序
	num=-1#判断是哪一种物质
	
	for j in range(len(path_list)):
		with open(path + path_list[j], 'r') as csvfile:
			x = np.loadtxt(csvfile)
			for i in range(len(x)):
				x[i]=float(x[i])

			x=np.mat(x)
			cos=cosine_similarity(x, y)
			#print(cos)
			if cos>0.9:
				num=j
				
	return num


def signal():
	# f=50 hz
	print('Generating simulated experiment')
	
	x=np.arange(0,1000,1)
	#x = np.arange(60, 1801, 1)
	g1=norm(loc = 100, scale = 1.0) # generate three gaussian as a signal
	g2=norm(loc = 200, scale = 3.0)
	g3=norm(loc = 750, scale = 5.0)
	g4 = norm(loc=550, scale=4.0)
	g5 = norm(loc=350, scale=8.0)
	
	signal=10000*(g1.pdf(x)+g2.pdf(x)+g3.pdf(x)+g4.pdf(x)+g5.pdf(x))
	baseline=10000*(0.015*np.sin(np.pi*x*5/x.max())+0.025*np.cos(np.pi*x*3/x.max()))+1000 # sinusoidal baseline
	noise=np.random.random(x.shape[0])/500
	
	y_origin=signal
	y=signal+baseline+noise

	# plt.plot(y)
	# plt.figure()
	# plt.plot(y_origin)
	# plt.show()
	return x,y_origin,y

def PLS(data):
    nmr=np.matrix(data)
    nmr=nmr.T
    nmr_sm10 = whittaker_smooth(nmr, 2E10, d=5)
    X=nmr-nmr_sm10
    for i in range(len(data)):
        if X[i]<0:
            X[i]=0

    # B=MMS(nmr)
    # for i in range(len(data)):
    #         nmr[i]=B[i]
    matplotlib.rcParams['font.sans-serif']=['SimHei']   # 用黑体显示中文
    matplotlib.rcParams['axes.unicode_minus']=False     # 正常显示负号
    return nmr,X,nmr_sm10

if __name__ == '__main__':
	import time
	
	start = time.clock()
	# long running
	# do something other
	
	
	
	# x,y_origin,fp=signal()
	# np.savetxt('./test/test1.txt', y_origin)
	# np.savetxt('./test/test.txt',fp)
	# z = baseline_als(10,0.00001,100,fp)
	# y2=fp-z
	# y2[y2 < 0] = 0
	# y2 = np.array(y2)
	# plt.plot(y2)
	# plt.show()
	# y22, matrix2 = fit(x,y2)
	# matrix2 = pd.DataFrame(matrix2)
	# matrix2.to_csv('./test/peak.csv', header=None)
	data_name='F:/laman/小拉曼/环己烷/c6-1s-1-50mw-8.csv'
	
	#data_name="F:/laman/生物大分子光谱数据/固体测试/DDT/raman/1μM/300mV-2s-1.csv"

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
	
	#y=np.loadtxt('F:/laman/test/高分子处理后数据/DDT-固体.txt')
	fp=data
	
	#fp=np.loadtxt('./test/test.txt')
	x,y=pre_data(xp,fp)
	np.savetxt('./test/test-nai.txt', y)
	# for i in range(335):  # 去除过小的杂乱信号
	# 	y[i] = 0
	yp = fp
	#y = np.mat(fp)
	#y = np.mat(fp)
	num = judge(y)
	#num = -1
	x = np.arange(65, 1701, 1)
	num=-2
	if num<0:
		n = [num]
		print('该物质的字典不在数据库中')
		np.savetxt('./test/num.txt', n)
		# A, y, C = PLS(yp)
		np.savetxt('./test/test.txt',yp)
		
		#yp=np.loadtxt('yuan.txt')
		#plt.plot(yp)
		#xx=np.arange(200, 3001, 1)
		yp = savgol_filter(yp,5, 3, mode="nearest")
		
		y2=baselineWavelet(yp)
		
		x,y2=pre_data(xp,y2)
		y2[y2 < 300] = 0
		# y111=np.zeros((1,500))
		# for i in range(500):
		# 	y111[0][i]=y2[i+1500]
		# y11, y22, noise = PLS(y111[0])
		# for i in range(500):
		# 	y2[i+1500]=y22[i]
		
		plt.plot(y2)
		plt.show()
		
		y2 = np.array(y2)
		y22, matrix2,p = fit(x,y2)
		plt.plot(y22)
		matrix2 = pd.DataFrame(matrix2)
		matrix2.to_csv('./test/peak.csv', header=None)
		plt.show()
	print(num)
	if num >= 0:
		n=[num]
		np.savetxt('./test/num.txt', n)
		print('该物质的字典在数据库中')
		
		x, yp = pre_data(xp, yp)
		np.savetxt('./test/test-nai.txt', yp)
		y2 = baselineWavelet(yp)
		y2[y2 < 0] = 0
		y2 = np.array(y2)
		

		y22, matrix2,p = fit(x,y2)
		
		path1 = "./标准峰/"  # 字典数据库
		path_list1 = os.listdir(path1)
		path_list1.sort()  # 对读取的路径进行排序
		with open(path1 + path_list1[num], 'r') as csvfile:
			peakp = np.loadtxt(csvfile)
		for i in range (len(p)):
			p[i]=p[i]+65
		peakp_real=np.zeros((1,len(peakp)))
		
		for j in range(len(peakp)):
			for i in range(len(p)):
				if abs(p[i]-peakp[j])<15:
					n=int(p[i]-65)
					intensity=y2[n]
					print(intensity)
					peakp_real[0][j]=intensity

		np.savetxt('./test/peakp.txt', peakp_real)
		path = "./lib/"  # 字典数据库
		path_list = os.listdir(path)
		path_list.sort()  # 对读取的路径进行排序
		df = pd.read_csv(path + path_list[num])
		data = pd.DataFrame(df)  # 将data的行列转换
		data.to_csv('./test/peak.csv', index=False,header=1)
	end = time.clock()
	print("代码执行时间是：")
	print(end - start)