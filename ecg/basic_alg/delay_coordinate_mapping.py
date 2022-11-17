"""
@Project ：
@File    ：delay_coordinate_mapping.py
@Author  ：DFQX
@Date    ：2022/11/15 22:45 
@Description: 使用delay-coordinate mapping算法
算法来源：A Real Time QRS Detection Using Delay-Coordinate Mapping for the
Microcontroller Implementation
"""

import utils.ecg_display as dp
import utils.signal_process as sig
import numpy as np
import data.read_data as rd


def inc_func(i, max_val):
    """
    自增函数
    """
    return i + 1 if i + 1 < max_val else 0


def dec_func(i, max_val):
    """
    自减函数
    """
    return i - 1 if i - 1 >= 0 else max_val - 1


def det_3(array):
    """
    计算3阶矩阵的行列式
    :param array: list(list), 方阵
    :return: 结果值
    """
    if len(array[0]) != 3 and 3 != len(array):
        raise '矩阵应该为3阶矩阵！'
    sum_arr, size = 0, len(array)
    for i in range(size):
        integral1, integral2 = 1.0, 1.0
        ai, aj, az = 0, i, i
        for j in range(size):
            integral1 *= array[ai][aj]
            integral2 *= array[ai][az]
            ai = inc_func(ai, size)  # 等价于 i + 1 if i + 1 < size else 0
            aj = inc_func(aj, size)
            az = dec_func(az, size)  # 等价于 i -1 if i-1 >= 0 else size-1
        sum_arr += integral1 - integral2
    return sum_arr


def det_2(array):
    """
    2阶矩阵行列式
    """
    if len(array) != 2 and len(array[0]) != 2:
        raise '矩阵应该为2阶矩阵！'
    return array[0][0] * array[1][1] - array[0][1] * array[1][0]


def max_func(arr):
    """
    一维列表中最大值和下标
    :param arr:  一维数组
    :return: 最大值， 最大值下标
    """
    m_v, m_i = 0, 0
    for idx, val in enumerate(arr):
        if val > m_v:
            m_v = val
            m_i = idx
    return m_v, m_i


def diff(data):
    result = []
    for idx in range(1, len(data)):
        result.append(data[idx] - data[idx - 1])
    return result


def delay_cor(ecg_data, delay_ms=9, fs=250):
    """
    延迟坐标算法
    :param delay_ms: int, 延时时长
    :param ecg_data: list, 心电数据
    :param fs: int, 采样率
    :return:
    """
    delay_n = int(delay_ms * fs / 1000 + 0.5)  # 延迟坐标点

    # ----------------滤波-------------------
    y_n = sig.bandpass_filter(ecg_data, fs, 1, 35)
    # -------------延迟坐标存储----------------
    x, y = [], []
    # -------------行列式存储----------------
    d = []  # d(n)矩阵
    v = []  # v(n)行列式值
    v_count = 0  # v计数
    # -----------阈值初始化(2秒)---------------
    thrs1, thrs2 = 0, 0  # eta_1 eta_2
    thrs1_arr, thrs2_arr = [], []
    # ------------llv和rlv存储----------------
    llv, rlv = [], []
    llv_q, rlv_q = [0 for val in range(5)], [0 for val in range(5)]
    q_i = 0  # llv_q 和 rlv_q的下标
    thrs_min_regions = []
    llv_sum, rlv_sum = [], []
    # --------------qrs波坐标----------------
    qrs, qrs_i = [], []  # r波波峰, r波坐标
    new_qrs = False
    refra_count, refra_tag = 1, False  # 不应期

    for i in range(delay_n, len(ecg_data)):
        x.append(ecg_data[i - delay_n])
        y.append(ecg_data[i])

        # -----------计算d(n)-----------
        idx = i - delay_n
        if idx >= 1:  # 最少有xi(0),xi(1)
            d.append([y[idx] - y[idx - 1], x[idx] - x[idx - 1]])
        if idx >= 2:  # len(d)>=2
            v.append(det_2(d[-2:]))  # np.linalg.det(d)
            v_count += 1

        # ------计算llv(n)和rlv(n)-------
        if idx >= 4:  # 最少需要3个v值
            llv.append(v[-2] + v[-3])  # llv和rlv会比v长度少2
            rlv.append(v[-2] + v[-1])

        # ---------阈值判断-------------
        if idx >= 6:
            llv_sum.append(llv[-1] + llv[-2] + llv[-3])  # llv_sum和rlv_sum会比v长度少4
            rlv_sum.append(rlv[-1] + rlv[-2] + rlv[-3])
            # ------使用两秒来进行初始化阈值----------
            if idx // fs < 2:
                thrs1 = 0.25 * max(llv_sum)
                thrs2 = thrs1
                if idx // fs < 1:   # 前一秒不进行阈值判断
                    continue

            if not new_qrs and not refra_tag and len(llv_sum) > 2:
                if abs(llv_sum[-1]) >= thrs1 and abs(rlv_sum[-1]) >= thrs1:  # case i
                    qrs.append(v[-1])
                    qrs_i.append(v_count - 1)
                    new_qrs = True
                    thrs1_arr.append(thrs1)
                    thrs2_arr.append(thrs2)
                elif abs(llv_sum[-1]) >= thrs2 or abs(rlv_sum[-1]) >= thrs2:  # case ii
                    qrs.append(v[-1])
                    qrs_i.append(v_count - 1)
                    new_qrs = True
                    thrs1_arr.append(thrs1)
                    thrs2_arr.append(thrs2)
                else:
                    new_qrs = False

        # --------------回找------------------
        if idx >= 6 and len(qrs_i) > 1 and not refra_tag:
            rr_i = qrs_i if len(qrs_i) <= 9 else qrs_i[-9:]  # 计算rr间期均值，一般8个rr间期的
            rr = diff(rr_i)
            rr_avg = sum(rr) / len(rr)
            if v_count - 1 - qrs_i[-1] > rr_avg * 1.5:  # 与上一次r波波峰位置大于最近的rr间期的1.5倍
                refra = int(0.2 * fs + 0.5)  # 不应期
                start = qrs_i[-1] + refra
                l_v, l_i = max_func(llv_sum[start:-1])
                r_v = max(rlv_sum[start:-1])
                r_i = start + l_i
                if l_v > thrs1 * 0.8 and r_v > thrs1 * 0.8:
                    qrs.append(v[r_i - 1])
                    qrs_i.append(r_i - 1)
                    new_qrs = True
                    thrs1_arr.append(thrs1 * 0.8)
                    thrs2_arr.append(thrs2 * 0.5)
                elif l_v > thrs2 * 0.5 or r_v > thrs2 * 0.5:
                    qrs.append(v[r_i - 1])
                    qrs_i.append(r_i - 1)
                    new_qrs = True
                    thrs1_arr.append(thrs1 * 0.8)
                    thrs2_arr.append(thrs2 * 0.5)
                else:
                    new_qrs = False

        # --------------更新阈值------------------
        if new_qrs and not refra_tag:
            llv_q[q_i] = llv_sum[-1]
            rlv_q[q_i] = rlv_sum[-1]
            q_i = inc_func(q_i, len(llv_q))
            thrs_min = min(min(llv_q), min(rlv_q))
            thrs_max = max(max(llv_q), max(rlv_q))
            thrs_min_region = (thrs_max + thrs_min) / 2
            thrs_min_regions.append(thrs_min_region)
            # 两个阈值最少需要两个或者三个thrs_min_region值
            if len(thrs_min_regions) >= 3:
                thrs1 = sum(thrs_min_regions[-2:]) / 2
                thrs2 = sum(thrs_min_regions[-3:]) / 3
            else:
                thrs1 = sum(thrs_min_regions) / (len(thrs_min_regions))
                thrs2 = thrs1
            refra_tag = True
            new_qrs = 0

        # -----------------不应期限制-------------
        if refra_tag:
            refra_count += 1
            if refra_count >= 0.2 * fs:  # 200ms不应期
                refra_tag = False
                refra_count = 1
    return v, llv_sum, rlv_sum, qrs, qrs_i, thrs1_arr, thrs2_arr


def d_cor(ecg_data, delay=9, fs=250):
    """
    坐标延迟，用作绘图研究
    """
    x_data, y_data = [], []
    d_point = int(delay * fs / 1000 + 0.5)
    x_data = ecg_data[0:-d_point]
    y_data = ecg_data[d_point:]
    dp.plot_ecg_delay(x_data, y_data)
    return x_data, y_data


def d_cor_3d(ecg_data, delay=9, fs=250):
    """
    坐标延迟，用作绘图研究
    """
    x_data, y_data, z_data = [], [], []
    d_point = int(delay * fs / 1000 + 0.5)
    x_data = ecg_data[0:-d_point * 2]
    y_data = ecg_data[d_point:-d_point]
    z_data = ecg_data[d_point * 2:]
    dp.plot_ecg_delay_3d(x_data, y_data, z_data)
    return x_data, y_data, z_data


def left_padding(data, num=4):
    """
    左填充0
    :param data: 数据
    :param num: 默认4
    :return: list
    """
    return [0 for v in range(4)] + data


if __name__ == '__main__':
    path = '../../data/mit-bih-arrhythmia'
    dat_name = '105.dat'
    # fs = 250
    fs = 360
    # ecg_data1, ecg_data2 = rd.read_f16(path, dat_name)
    ecg_data1, ecg_data2 = rd.read_f212(path, dat_name)
    ecg = ecg_data1[:100000]
    # ecg1 = sig.bandpass_filter(ecg, fs, 1, 35)
    # v_n = delay_cor(ecg1, delay_ms=9, fs=fs)
    # d_cor(ecg, 9, fs)  # 绘制坐标延迟
    # d_cor_3d(ecg1, 9, fs)  # 绘制3D坐标
    v, llv_sum, rlv_sum, qrs, qrs_i, thrs1_arr, thrs2_arr = delay_cor(ecg, 9, fs)

    dp.plot_simple_comp(ecg, v)
    # dp.plot_simple_comp(ecg, ecg1)
    # dp.plot_simple_comp(SUM_LLV_LIST, SUM_RLV_LIST)
    # dp.plot_simple_comp(ecg3, ecg4)
    dp.plot_peak_dot(v, qrs_i, qrs)
    # dp.plot_peak_dot_th1_th2(V, QRS_i, QRS, thrs1_arr, thrs2_arr)

    dp.plot_peak_dot_llv_rlv(v, qrs_i, qrs, left_padding(llv_sum), left_padding(rlv_sum), thrs1_arr, thrs2_arr)
