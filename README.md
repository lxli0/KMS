# K–M slope (KMS)

K–M slope (KMS) is a parameter derived from Visibility Graph Analysis (Lacasa et al., 2008 PNAS; Telesca et al., 2013 Physica A).
We verify that there is a positive between KMS and b-value from Gutenberg–Richter law. Thus, it can be used to estimate the b-value.

Here, We provide:
1.  The Matlab code to calculate KMS_0 (KMS/b ratio) under a given catalog size. --KMS2b_ratio.m
2.  The Matlab code to calculate the temporal evolution of the b-value by KMS and the improved KMS methods. The earthquake data in the 
    northeastern Tibetan Plateau from 1980 to 2020 is the input of the code and can be freely replaced. --Temporal_bvalue_by_KMS.m
3.  A relatively fast (compared with the above two scrips, the time complexity is O(n^2)) Matlab code for calculating the b-value
    of a given magnitude sequence based on the improved KMS method. The script is still based on the classic definition of a
    visibility graph. --KMS_Improved_tradditional.m
5.  An efficient (the time complexity is O(n*logn)) Matlab code for calculating the b-value of a given magnitude sequence based on
    the improved KMS method. The script is based on the Divide & Conquer (DC) strategy (Lan et al., 2015). --KMS_Improved_DC_algorithm.m

Reference: 
1.  Linxuan Li, Gang Luo, and Mian Liu (2023). The K–M slope: a potential supplement for b-value. Seismological Research Letters 94, 1892-1899. doi: 10.1785/0220220268
2.  Linxuan Li and Gang Luo (submitted). Can we obtain reliable seismic b-values for real-time catalogs?
3.  Lucas Lacasa, Bartolo Luque, Fernando Ballesteros, Jordi Luque, and Juan Carlos Nuño (2008). From time series to complex networks: The visibility graph. Proceedings of the National Academy of Sciences 105(13), 4972-4975. doi: 10.1073/pnas.0709247105
4.  Xin Lan, Hongming Mo, Shiyu Chen, Qi Liu, and Yong Deng (2015). Fast transformation from time series to visibility graphs. Chaos: An Interdisciplinary Journal of Nonlinear Science. 25:083105. doi: 10.1063/1.4927835


Contact: Linxuan Li.  lli7@caltech.edu or lucas.linxuan.li@gmail.com
