import pandas
import sys
import numpy as np

# 从命令行参数获取输入和输出文件的路径
infile = sys.argv[1]
recfile = sys.argv[2]

# 读取输入文件，假设文件以制表符分隔且没有列标题
inf = pandas.read_csv(infile, sep='\t', header=None)
# 设置DataFrame的索引为第二列（索引为1）
inf.index = inf[1]

# 创建一个新的DataFrame来存储重组率数据
rec = pandas.DataFrame()
rec['start.pos'] = list(inf[1])
rec['recom.rate.perbp'] = 0.0
rec['recom.rate.cMperMb'] = 0.0
rec['gm'] = list(inf[2])
rec['next_gm'] = rec['gm'].shift(periods=-1)
# 设置最后一行的next_gm为其gm值，避免NaN
rec.iloc[-1, -1] = rec.iloc[-1, -2]
rec['next_pos'] = rec['start.pos'].shift(periods=-1)
# 将最后一个next_pos设置为0，避免除零错误
rec.iloc[-1, -1] = 0.0
# 计算重组率
rec['recom.rate.perbp'] = ((rec['next_gm'] - rec['gm']) / (rec['next_pos'] - rec['start.pos'])) / 100.0       #Mperbp
rec['recom.rate.cMperMb'] = (rec['next_gm'] - rec['gm']) / ((rec['next_pos'] - rec['start.pos']) / 1000000 )  #cMperMb
# 删除不再需要的列
rec.drop(['next_gm', 'next_pos'], axis=1, inplace=True)
# 设置最后一个recom.rate.perbp和recom.rate.cMperMb为0
rec.iloc[-1, -2] = 0.0
rec.iloc[-1, -3] = 0.0

# 将结果输出到文件
rec.to_csv(recfile, sep=' ', index=None, float_format='%.15f')

