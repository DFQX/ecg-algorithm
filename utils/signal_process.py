"""
@Project ：
@File    ：signal_process.py
@Author  ：DFQX
@Date    ：2022/11/15 22:55 
@Description: 信号处理相关方法
"""
from scipy import signal


def bandpass_filter(data, fs, low=5, high=15):
    """
    带通滤波
    :param data: list, 输入信号数据
    :param fs: int, 采样率
    :param low: int, 截止频率1
    :param high: int, 截止频率2
    :return: list, 滤波后信号
    """
    low = 2 * low / fs
    high = 2 * high / fs
    b, a = signal.butter(3, [low, high], 'bandpass')
    filter = signal.filtfilt(b, a, data)  # data为要过滤的信号
    return filter
