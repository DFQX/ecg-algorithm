"""
@Project ：
@File    ：ecg_display.py
@Author  ：DFQX
@Date    ：2022/9/18 13:09 
@Description: 使用matplotlib打印心电图
"""

from matplotlib import pyplot as plt
import data.read_data as rd
import data.covert as cov

'''
函数：plt_ecg
描述：用于绘制只有一条导联数据的心电图，并且时长不能太长
输入：data:心电数据；fs:采样率；gain:增益; path:保存图片路径，默认不保存
输出：一张心电图，没有保存为文件，如果需要，使用plt.savefig()
'''


def plt_ecg(data, fs=250, gain=1, path=''):  # 打印ECG图纸,一个心电图
    """
    用于绘制只有一条导联数据的心电图，并且时长不能太长
    :param data: list, 心电数据
    :param fs: int, 采样率
    :param gain: int, 增益
    :param path: str， 默认为空，不用保存，设置路径后保存图片
    :return: 无
    """
    x = [i / fs for i in range(0, len(data))]       # 计算x轴时间
    y = [val / gain for val in data]                # 根据增益计算幅值
    plt.figure(figsize=(len(data)*25//(fs*14), 4))  # 画板大小固定，不要更改
    plt.xlabel("Time: s", fontsize=14)
    plt.ylabel("Voltage: mV", fontsize=14)
    plt.margins(x=0)
    ax = plt.gca()
    # 心电图纸：横坐标一般小格0.04s,大格0.2s; 纵坐标小格0.1mv，大格0.5mv
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))   # 设置x轴主刻度
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.04))  # 设置x轴次刻度
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))   # 设置y轴主刻度
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))   # 设置y轴次刻度
    # 绘制大格和小格的线条
    ax.grid(which='major', axis="both", linewidth=0.75, linestyle='-', color='r')
    ax.grid(which='minor', axis="both", linewidth=0.25, linestyle='-', color='r')
    plt.ylim([-2.5, 2.5])   # y轴值为-2.5~2.5mV之间
    plt.plot(x, y, 'black', linewidth=0.9)   # 心电图形绘制
    plt.savefig(path) if path.strip() != '' else None
    plt.show()


if __name__ == '__main__':
    path = '../data/paf'
    dat_name = 'n01.dat'
    ecg_data1, ecg_data2 = rd.read_f16(path, dat_name)
    ecg_data1 = cov.covert_freq(ori_freq=250, tar_freq=200, data=ecg_data1)
    plt_ecg(ecg_data1[:1001], fs=200, gain=200, path='1.jpg')

if __name__=='__main1__':
    plt.figure(figsize=(10,4))  # 画板大小固定，不要更改
    ax = plt.gca()
    # 心电图纸：横坐标一般小格0.04s,大格0.2s; 纵坐标小格0.1mv，大格0.5mv
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))  # 设置x轴主刻度
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.04))  # 设置x轴次刻度
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))  # 设置y轴主刻度
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))  # 设置y轴次刻度
    # 绘制大格和小格的线条
    ax.grid(which='major', axis="both", linewidth=0.75, linestyle='-', color='r')
    ax.grid(which='minor', axis="both", linewidth=0.25, linestyle='-', color='r')
    plt.ylim([-2.5, 2.5])  # 纵坐标范围
    plt.xlim([0, 5])       # 横坐标范围
    plt.show()
