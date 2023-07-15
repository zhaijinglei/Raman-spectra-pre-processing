

import numpy as np
import math
import matplotlib.pyplot as plt

from baselineCorrection import fit

def snr_out(signal_source, signal_source_noise):
    signal_noise = signal_source - signal_source_noise
    mean_signal_source = np.mean(signal_source)
    signal_source = signal_source - mean_signal_source
    snr = 10 * math.log(np.sum(signal_source ** 2) / np.sum(signal_noise ** 2), 10)
    return snr
    
    return (snr)


def MMS(data):  # 最大值归一化
    return [float(i) / (max(data)) for i in data]


def get_rmse(records_real, records_predict):
    """
    均方误差 估计值与真值 偏差
    """
    if len(records_real) == len(records_predict):
        return math.sqrt(sum([(x - y) ** 2 for x, y in zip(records_real, records_predict)]) / len(records_real))
    else:
        return None


def cosine_similarity(x, y):  # 余弦相似度
    num = x.dot(y.T)
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    n = num / denom
    if n > 1: n = 1.0
    return n



SN = np.loadtxt('./test/signal.txt')
S1 = np.loadtxt('./test/omp_signal1.txt')#静态字典
S2 = np.loadtxt('./test/omp_signal3.txt')#动态字典
S3 = np.loadtxt('./test/ pf.txt')
S4 = np.loadtxt('./test/ann.txt')
yuan = np.loadtxt('./test/ y_origin.txt')
plt.plot(yuan)
plt.figure()
plt.plot(S1)
plt.figure()
plt.plot(S2)
plt.show()

x = np.arange(70,1101,1)
S3=fit(x,SN,3)
np.savetxt('./test/ pf.txt',S3)
SNR1 = snr_out(SN, yuan)#信噪比，越大越好
SNR2 = snr_out(S1, yuan)
SNR3 = snr_out(S2, yuan)
SNR4= snr_out(S3, yuan)
SNR5= snr_out(S4, yuan)

rmse1 = get_rmse(SN, yuan)#均方根误差，越小越好
rmse2 = get_rmse(S1, yuan)
rmse3 = get_rmse(S2, yuan)
rmse4 = get_rmse(S3, yuan)
rmse5 = get_rmse(S4, yuan)

cos1 = cosine_similarity(SN, yuan)#余弦相似度，越小越好
cos2 = cosine_similarity(S1, yuan)
cos3 = cosine_similarity(S2, yuan)
cos4 = cosine_similarity(S3, yuan)
cos5 = cosine_similarity(S4, yuan)

# plt.plot(yuan, '-k')
# plt.figure()
# plt.plot(SN, '-r')
# plt.figure()
# plt.plot(S, '-b')
# plt.figure()
# plt.plot((S-yuan),'-k')
# plt.show()
from scipy.stats import pearsonr

pccs1 = pearsonr(SN, yuan)#皮尔逊相关系数
pccs2 = pearsonr(S1, yuan)
pccs3 = pearsonr(S2, yuan)
pccs4 = pearsonr(S3, yuan)
pccs5 = pearsonr(S4, yuan)

print(SNR1, SNR2,SNR3,SNR4,SNR5)
print(rmse1, rmse2,rmse3,rmse4,rmse5)
print(cos1, cos2,cos3,cos4,cos5)
print(pccs1[0], pccs2[0],pccs3[0],pccs4[0],pccs5[0])