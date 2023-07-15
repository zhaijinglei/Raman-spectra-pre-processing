import numpy as np
from lmfit.models import (VoigtModel)  # pip install lmfit
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os
from pandas.core.frame import DataFrame
from scipy import interpolate
from make_dictionary import fit
from baselinewavelet import baselineWavelet
from scipy.sparse import csc_matrix, eye, diags
# from baselineCorrection import baseline_als
from scipy.stats import norm


def mean(data):
	""""
	data 横轴方向表示重复测量的数据，纵轴方向表示波段数量
	"""
	# 计算平均光谱，实际就是x值
	s_mean = np.mean(data, axis=1)
	
	for i in range(len(s_mean)):
		if s_mean[i] < 0:
			s_mean[i] = 0
	return s_mean


def read():
	path = "F:/laman/生物大分子光谱数据/液体sers/DDT/1μM/"  # 待读取的文件夹
	path_list = os.listdir(path)
	path_list.sort()  # 对读取的路径进行排序
	with open(path + path_list[0], 'r') as csvfile:
		reader = csv.reader(csvfile)
		column = [row[0] for row in reader]
	
	with open('./高分子数据/information.csv', 'a', encoding='utf-8', newline='') as fp:
		# 写
		writer = csv.writer(fp)
		writer.writerow(column)
	
	for i in range(20):
		csv_path = path + path_list[i]
		with open(csv_path, 'r') as csvfile:
			reader = csv.reader(csvfile)
			column = [row[1] for row in reader]
		
		with open('./高分子数据/information.csv', 'a', encoding='utf-8', newline='') as fp:
			# 写
			writer = csv.writer(fp)
			writer.writerow(column)
	df = pd.read_csv('./高分子数据/information.csv')
	data = df.values  # data是数组，直接从文件读出来的数据格式是数组
	index1 = list(df.keys())  # 获取原有csv文件的标题，并形成列表
	data = list(map(list, zip(*data)))  # map()可以单独列出列表，将数组转换成列表
	data = pd.DataFrame(data, index=index1)  # 将data的行列转换
	data.to_csv('./高分子数据/information1.csv', header=1)
	data.to_csv('./高分子数据/information5.csv', header=1)


def wjm():
	read()
	with open("F:/laman/生物大分子光谱数据/液体sers/DDT/1μM/DATA-215154-X0-Y1-2590.csv", 'r') as csvfile:
		reader = csv.reader(csvfile)
		
		x1 = [row[0] for row in reader]


	txt_path = r"./高分子数据/information1.csv"
	df = pd.read_csv(txt_path, header=0, index_col=0, engine="python", encoding="utf-8")
	# 行列名称
	df_arr = df.values
	s_mean = mean(df_arr)
	for i in range(len(s_mean)):
		if s_mean[i] < 0:
			s_mean[i] = 0

	with open('./高分子数据/information2.csv', 'a', encoding='utf-8', newline='') as fp:
		# 写
		writer = csv.writer(fp)
		writer.writerow(x1)
		writer.writerow(s_mean)
	df = pd.read_csv('./高分子数据/information2.csv')
	data = df.values  # data是数组，直接从文件读出来的数据格式是数组
	index1 = list(df.keys())  # 获取原有csv文件的标题，并形成列表
	data = list(map(list, zip(*data)))  # map()可以单独列出列表，将数组转换成列表
	data = pd.DataFrame(data, index=index1)  # 将data的行列转换
	data.to_csv('./高分子数据/information3.csv', header=1)
	txt_path = r"./高分子数据/information3.csv"
	df = pd.read_csv(txt_path, header=0, index_col=0, engine="python", encoding="utf-8")
	# 行列名称
	columns_n = df.columns.values.tolist()
	row_n = df._stat_axis.values.tolist()
	df_arr = df.values
	dabai1_arr_data_you1 = DataFrame(df_arr, index=row_n, columns=columns_n)
	os.remove('./高分子数据/information.csv')
	os.remove('./高分子数据/information1.csv')
	os.remove('./高分子数据/information2.csv')
	os.remove('./高分子数据/information3.csv')
	return dabai1_arr_data_you1, x1, s_mean
def pre_data(x, y):
	f1 = interpolate.interp1d(x, y, kind='linear')  # 插值模块
	x_pred = np.linspace(200, 2200, num=2001)
	y_pred = f1(x_pred)

	return x_pred, y_pred

dabai1_arr_data_you1, x2, smean=wjm()
for i in range(len(x2)):
	x2[i]=float(x2[i])
	smean[i] = float(smean[i])
x,y=pre_data(x2,smean)
np.savetxt('./高分子处理后数据/DDT-液体.txt',y)