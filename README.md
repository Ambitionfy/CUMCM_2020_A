# cumcm_2020_A
2020全国大学生数学建模大赛A题
原推国家级一等奖，因附录格式问题，降为省级一等奖

本题为回焊炉问题
1.通过建立基于热传导与热对流的模型，可以很好的模拟回焊炉中芯片的温度变化。
2.对于热传导形成的三对角矩阵，采用LU分解，通过不断迭代来求解每时刻的温度变化情况。
3.对于参数的选择，采用进化算法（这里采用粒子群），来精确求解各参数。

本题遗憾：对题意理解仍稍有偏颇
本题中有句：回焊炉开启后，在极短时间内达到稳定，我们简单的理解为在很短时间内的温度变化小于一定阈值，即认为稳定。而没有更直接的理解为温区呈线性。
因此，本模型中求出的回焊炉间隙中，温区温度呈高元曲线，芯片在该区域会稍有温度下降。
