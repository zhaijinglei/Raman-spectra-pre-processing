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


def pre_data(x, y):
	f1 = interpolate.interp1d(x, y, kind='linear')  # 插值模块
	x_pred = np.linspace(70, 1100, num=1031)
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


def signal_genarate():
	# f=50 hz
	print('Generating simulated experiment')
	data_name = 'F:/laman/小拉曼/硫/硫-1s-1-150mw-10.csv'
	data=np.loadtxt('./硫.txt')
	
	# with open(data_name, 'r') as csvfile:
	# 	reader = csv.reader(csvfile)
	# 	data = [row[1] for row in reader]  #
	# for i in range(len(data)):
	# 	data[i] = float(data[i])
	with open(data_name, 'r') as csvfile:
		reader = csv.reader(csvfile)
		xp = [row[0] for row in reader]  #
	for i in range(len(xp)):
		xp[i] = float(xp[i])
		
	df = pd.read_csv('./data_noise.csv', header=None)
	data_noise = df.values
	noise=np.zeros((1,len(data_noise[0])))
	
	for i in range(len(data_noise)):
		noise=noise+50*data_noise[i]
	signal = data
	# x, y = pre_data(xp, signal)
	y=baselineWavelet(signal)
	x=np.arange(70,1101,1)
	baseline=100000*(0.015*np.sin(np.pi*x*5/x.max())+0.025*np.cos(np.pi*x*3/x.max()))+1000+3000 # sinusoidal baseline
	np.savetxt('./test/background.txt',baseline)
	plt.figure()
	plt.plot(baseline)
	plt.figure()
	y_origin=y
	signal=y+baseline+noise[0]
	
	plt.plot(x, y_origin)
	plt.figure()
	plt.plot(x,signal)
	plt.show()
	return x,y_origin,signal

if __name__ == '__main__':
	

	x,y_origin,signal=signal_genarate()
	
	np.savetxt('./test/ y_origin.txt', y_origin)
	np.savetxt('./test/signal.txt',signal)
	
	# z = baseline_als(10,0.00001,100,signal)
	# y2=signal-z
	y2=baselineWavelet(signal)
	y2[y2 < 0] = 0

	
	y2 = np.array(y2)
	plt.plot(x, y2)
	plt.figure()
	y22, matrix2,p = fit(x,y2)
	bizhi=(max(y2))/(max(y22))
	print(bizhi)
	plt.plot(x,y22)
	plt.show()
	matrix2 = pd.DataFrame(matrix2)
	matrix2.to_csv('./test/peak_sti.csv', header=None)
	