#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：ecg-algorithm
@File    ：read_data.py
@Author  ：DFQX
@Date    ：2022/8/10 17:59
@Description: 读取ECG数据
"""
import os
import matplotlib.pyplot as plt


def read_hea(atr_path, hea_name):
    """
    读取注释文件'.hea'的内容
    :param atr_path: str, 注释文件路径
    :param hea_name: str, 注释文件名
    :return: dict
    """
    data = {}
    with open(os.path.join(atr_path, hea_name), 'r') as file:
        line1 = file.readline().strip().split(' ')
        data['filename'] = line1[0]
        data['lead_num'] = int(line1[1])
        data['sample_freq'] = int(line1[2])
        data['sample_num'] = int(line1[3])
        # 分别是信号格式，增益，采样精度，零值，第一个值（用于偏差校验）
        sformat, gain, bit_res, zero_value, first_value = [0] * data['lead_num'] \
            , [0] * data['lead_num'], [0] * data['lead_num'] \
            , [0] * data['lead_num'], [0] * data['lead_num']
        for lead_idx in range(data['lead_num']):  # 每个导联进行读取
            line = file.readline().strip().split(' ')
            sformat[lead_idx], gain[lead_idx], bit_res[lead_idx] = int(line[1]) \
                , int(line[2]) if '/' not in line[2] else int(line[2][:-3]), int(line[3])
            zero_value[lead_idx], first_value[lead_idx] = int(line[4]), int(line[5])
        data['sformat'] = sformat
        data['gain'] = gain
        data['bit_res'] = bit_res
        data['zero_value'] = zero_value
        data['first_value'] = first_value
        return data


def read_f212(file_path, file_name):
    """
    读取 format 212 格式文件 '.dat'
    :param file_path: str, 数据路径
    :param file_name: str, 数据文件名
    :return:
    """
    data = []
    with open(os.path.join(file_path, file_name), 'rb') as file:
        data = file.read()  # 读取二进制文件
    ecg_data_lead1, ecg_data_lead2 = [], []
    for i in range(0, len(data), 3):
        lead1 = ((data[i + 1] & 0x0f) << 8) + (data[i])  # 导联1数据读取
        lead2 = ((data[i + 1] & 0xf0) << 4) + (data[i + 2])  # 导联2数据读取
        if lead1 > 2047:  # 由于加载的数据是无符号型，转换为二进制补码形式：value>2^11-1 为负值
            lead1 -= 4096
        if lead2 > 2047:
            lead2 -= 4096
        ecg_data_lead1.append(lead1)
        ecg_data_lead2.append(lead2)
    # plot_data(ecg_data_lead2)  # 打印导联1数据
    return ecg_data_lead1, ecg_data_lead2


def read_f16(file_path, file_name):
    """
    读取format 16格式的文件
    :param file_path: str, 文件路径
    :param file_name: str, 文件名
    :return: list, list, 导联1和导联2的数据
    """
    ecg_data = []
    with open(os.path.join(file_path, file_name), 'rb') as file:
        ecg_data = file.read()
    ecg_data_lead1, ecg_data_lead2 = [], []
    for idx in range(0, len(ecg_data), 4):
        lead1 = (ecg_data[idx + 1] << 8) + ecg_data[idx]
        lead2 = (ecg_data[idx + 3] << 8) + ecg_data[idx + 2]
        if lead1 > 2 ** 15 - 1:
            lead1 -= 2 ** 16
        if lead2 > 2 ** 15 - 1:
            lead2 -= 2 ** 16
        ecg_data_lead1.append(lead1)
        ecg_data_lead2.append(lead2)
    # plot_data(ecg_data_lead1)  # 绘制导联1的心电图
    # plot_data(ecg_data_lead2)  # 绘制导联2的心电图
    return ecg_data_lead1, ecg_data_lead2


def read_atr(atr_path, atr_name):
    """
    读取注释文件
    :param atr_path: str, 文件路径
    :param atr_name: str, 文件名
    :return: list, list, 注释， 时间间隔
    annotation的相关信息参见：https://www.physionet.org/physiotools/wag/annot-5.htm
    """
    data = []
    with open(os.path.join(atr_path, atr_name), 'rb') as file:
        data = file.read()
    atr_time, annot = [], []
    atr_info, atr_info1, atr_info2 = [], [], []  # 存储两行按照字节读取的信息
    for idx, v in enumerate(data):  # 数据按照[[1,3,5,...],[2,4,6,...]]进行存储
        if idx % 2 == 0:
            atr_info1.append(v)
        else:
            atr_info2.append(v)
    atr_info.append(atr_info1)  # 存到一个列表中
    atr_info.append(atr_info2)
    idx = 0
    while idx < len(atr_info[0]):  # 根据MIT-BIH format来解析数据
        annoth = atr_info[1][idx] >> 2
        if annoth == 59:
            annot.append(atr_info[1][idx + 3] >> 2)
            atr_time.append(atr_info[0][idx + 2]
                            + (atr_info[1][idx + 2] << 8)
                            + (atr_info[0][idx + 1] << 16)
                            + (atr_info[1][idx + 1] << 24))
            idx += 3
        elif annoth == 60:
            pass
        elif annoth == 61:
            pass
        elif annoth == 62:
            pass
        elif annoth == 63:
            hilfe = ((atr_info[1][idx] & 3) << 8) + atr_info[0][idx]
            hilfe = hilfe + hilfe % 2
            idx += hilfe // 2
        else:
            atr_time.append(((atr_info[1][idx] & 0x3) << 8) + atr_info[0][idx])  # 低10位为时间间隔
            annot.append(atr_info[1][idx] >> 2)  # 高6位为注释
        idx += 1
    return annot, atr_time

def read_qrs(qrs_path, qrs_name):
    """
    读取.qrs文件，比如PAF数据库中，xxx.qrs文件
    :param qrs_path: str, 文件路径
    :param qrs_name: str, 文件名
    :return:
    """
    data = []
    with open(os.path.join(qrs_path, qrs_name), 'rb') as file:
        data = file.read()
    print(data[:8])


def plot_data(data):
    plt.figure(figsize=(12, 4))
    plt.plot(data)
    plt.show()


def plot_ecg_and_annot(data, annot, atr_time):
    """
    绘制心电信号和注释
    :param data: list, 心电信号
    :param annot: list, 注释
    :param atr_time: list, 注释和上一个注释之间的间隔
    :return: 无
    """
    plt.figure(figsize=(12, 4))
    plt.plot(data)
    for idx, v in enumerate(atr_time):
        if idx == 0:
            continue
        atr_time[idx] = atr_time[idx - 1] + atr_time[idx]
    r_point = [data[v] for v in atr_time]
    plt.plot(atr_time, r_point, 'r*')
    plt.show()


if __name__ == '__main__':
    path = './paf'
    dat_name = 'n01c.dat'
    hea_name = 'n01c.hea'
    atr_name = 'n01c.qrs'
    atr_info = read_hea(path, hea_name)
    ecg_data1, ecg_data2 = read_f16(path, dat_name)
    # 第一值校验
    if ecg_data1[0] != atr_info['first_value'][0] or ecg_data2[0] != atr_info['first_value'][1]:
        raise Exception('inconsistency in the first bit values')
    # # 减0值后除以增益，转换成mV单位
    ecg_data1 = [(v - atr_info['zero_value'][0]) / atr_info['gain'][0] for v in ecg_data1]
    ecg_data2 = [(v - atr_info['zero_value'][1]) / atr_info['gain'][1] for v in ecg_data2]
    # annot, atr_time = read_qrs(path, atr_name)
    annot, atr_time = read_atr(path, atr_name)
    plot_data(ecg_data1)
    # plot_ecg_and_annot(ecg_data1[5:], annot, atr_time)    # 打印ecg图和R波位置


# if __name__ == '__main__':
#     path = './mit-bih-arrhythmia'
#     dat_name = '100.dat'
#     hea_name = '100.hea'
#     atr_name = '100.atr'
#     atr_info = read_hea(path, hea_name)
#     ecg_data1, ecg_data2 = read_f212(path, dat_name)
#     # 第一值校验
#     if ecg_data1[0] != atr_info['first_value'][0] or ecg_data2[1] != atr_info['first_value'][1]:
#         raise Exception('inconsistency in the first bit values')
#     # 减0值后除以增益，转换成mV单位
#     ecg_data1 = [(v - atr_info['zero_value'][0]) / atr_info['gain'][0] for v in ecg_data1]
#     ecg_data2 = [(v - atr_info['zero_value'][1]) / atr_info['gain'][1] for v in ecg_data2]
#     annot, atr_time = read_atr(path, atr_name)
#     print(annot[:15])
#     plot_ecg_and_annot(ecg_data1, annot, atr_time)
