"""
@Project ：
@File    ：wfdb_read.py
@Author  ：DFQX
@Date    ：2022/9/17 16:06 
@Description:  使用wfdb进行读取ecg信号
"""
import wfdb
# './paf/n01'是下载来的PAF数据库的文件, channels是选择的导联
record = wfdb.rdrecord('./aha/0201')
annotation = wfdb.rdann('./aha/0201', 'atr')
wfdb.plot_wfdb(record=record, annotation=annotation, title='PAF n01', time_units='seconds')

