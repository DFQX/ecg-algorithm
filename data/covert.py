"""
@Project ：
@File    ：covert.py
@Author  ：DFQX
@Date    ：2022/9/18 14:50 
@Description:  一些转换工具
"""
import numpy as np


def covert_freq(ori_freq, tar_freq, data):
    """
    将信号频率进行转换, 只能进行降采样，这里不做插值
    即tar_freq<=ori_freq
    :param ori_freq: int, 原始信号频率
    :param tar_freq: int, 转换后信号频率
    :param data: list, 信号数据
    :return: list, 返回后的信号
    """
    if len(data) == 0:
        raise Exception('数据长度为0')
    if ori_freq < tar_freq:
        raise Exception('原始频率应该大于目标频率！')
    ori_time = np.arange(0, len(data), 1 / ori_freq)  # 原始数据坐标时间
    tar_time = np.arange(0, int(len(data) * tar_freq // ori_freq), 1 / tar_freq)  # 转换数据坐标时间
    tar_data, j = [], 0
    for idx, v in enumerate(ori_time):
        if idx >= len(data):
            break
        if tar_time[j] <= ori_time[idx]:
            tar_data.append(data[idx])
            j += 1
    return tar_data
