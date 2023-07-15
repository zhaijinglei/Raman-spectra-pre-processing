import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from numpy import *
import peakutils
from scipy.sparse.linalg import spsolve
from scipy import sparse
import pandas as pd

bW= importr('baselineWavelet')

def baselineWavelet(Spectrum_data):
    S=Spectrum_data
    # Transform the original data into R language format
    rFS=ro.FloatVector(S)
    # Creating a R array correspond to the data's dimension
    sS=ro.FloatVector(range(1,shape(S)[0]+1))
    ## Start the wavelet analysis from R's baselineWavelet package
    scales=ro.r.seq(1,127,1)
    cwtargs={'scales':scales,'wavelet':'mexh'}
    
    # Wavelet coeffecients calculated by cwt
    ####******************########
    wCoefs=ro.r.cwt(rFS,**cwtargs)
    ####******************########
    
    # Get the ridge of wavelet coeffecients and output it into png file.
    localMax=bW.getLocalMaximumCWT(wCoefs)
    ####******************########
    ridgeList=bW.getRidge(localMax, gapTh=2, skip=2)
    # Calculating the background contribution
    MPargs={'SNR.Th':1,'ridgeLength':5}
    # Peak infomation
    majorPeakInfo = bW.identifyMajorPeaks(rFS, ridgeList, wCoefs,scales,**MPargs)
    peakWidth=bW.widthEstimationCWT(rFS,majorPeakInfo)
    
    bCargs={'lambda':10,'differences':1}
    ####******************########
    backgr = bW.baselineCorrectionCWT(rFS,peakWidth,**bCargs)
    corrected=ro.FloatVector(array(rFS)-array(backgr))
    ####******************########
    #wavelet coeffecients, baseline, baseline-corrected spectrumarray(wCoefs),array(backgr),
    return array(corrected)

if __name__=='__main__':
    import csv
    # df=pd.read_csv('F:/laman/print/displaydiagram/c1/data2_2/information.csv')
    # base= df.values  # data是数组，直接从文件读出来的数据格式是数组


    with open('./test/test.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        fp = [row[0] for row in reader]  #

    for i in range(len(fp)):
        fp[i] = float(fp[i])

    y2 = baselineWavelet(fp)
    # base = pd.DataFrame(base)  # 将data的行列转换
    # print(base[0])
    # base[0]=baselineWavelet(base[0])
    
    # for i in range(10):
    #     z0 = baseline_als(10, 0.01, 10, base[i])
    #     base[i] = base[i] - z0
    #     base[i][base[i] < 0] = 0
    # np.set_printoptions(threshold=np.inf)
    # base = pd.DataFrame(base)  # 将data的行列转换
    # base.to_csv('F:/laman/print/displaydiagram/最小二乘校正.csv', header=0)
    # x = np.arange(0, 2048, 1)
    # plt.title('智能拉曼光谱荧光背景扣除算法')
    # plt.plot(x, base, '-r', label="校正后信号")
    # plt.plot(x, A, '-k', label="原始信号")
    # plt.legend(loc='upper left')
    # matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 用黑体显示中文
    # matplotlib.rcParams['axes.unicode_minus'] = False  # 正常显示负号
    # plt.show()

