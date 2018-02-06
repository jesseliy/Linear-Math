# TODO 实现 Gaussain Jordan 方法求解 Ax = b

"""
Draft submitted in 20180206.
Debug etc. NEEDED.
Revision of other def in  Jupyter Notebook are also needed
"""

""" Gaussian Jordan 方法求解 Ax = b.
    参数
        A: 方阵
        b: 列向量
        decPts: 四舍五入位数，默认为4
        epsilon: 判读是否为0的阈值，默认 1.0e-16

    返回列向量 x 使得 Ax = b
    返回None，如果 A，b 高度不同
    返回None，如果 A 为奇异矩阵
"""

"""
步骤1 检查A，b是否行数相同

步骤2 构造增广矩阵Ab

步骤3 逐列转换Ab为化简行阶梯形矩阵 中文维基链接

对于Ab的每一列（最后一列除外）
    当前列为列c
    寻找列c中 对角线以及对角线以下所有元素（行 c~N）的绝对值的最大值
    如果绝对值最大值为0
        那么A为奇异矩阵，返回None (你可以在选做问题2.4中证明为什么这里A一定是奇异矩阵)
    否则
        使用第一个行变换，将绝对值最大值所在行交换到对角线元素所在行（行c） 
        使用第二个行变换，将列c的对角线元素缩放为1
        多次使用第三个行变换，将列c的其他元素消为0

步骤4 返回Ab的最后一列
"""
def matxRound(M, decPts=4):
    for i in range(len(M)):
        interv = [round(x,decPts) for x in M[i]]
        M[i] = interv
	return M
	
def swapRows(M, r1, r2):
    M[r1],M[r2] = M[r2],M[r1]
	return M


def scaleRow(M, r, scale):
    if scale == 0:
        raise ValueError
    else:
        M[r] = [x*scale for x in M[r]]
	return M


def addScaledRow(M, r1, r2, scale):
    M[r1] = [x+y*scale for x,y in zip(M[r1],M[r2])]
	return M


def augmentMatrix(A, b):
    N = [x+y for x,y in zip(A,b)]
    return N

def gj_Solve(A, b, decPts=4, epsilon = 1.0e-16):

    if len(A[0]) != len(b):  # 如果 A，b 高度不同,返回None
        return None
    Ab = augmentMatrix(A,b)  # 增广矩阵
    col = len(Ab)
	row = len(Ab[0])
	
    for r in range(row):
		ind_max = r
		num_max = abs(Ab[r][r])
		for c in range(r,col):       # 寻找对角线及以下所有元素的绝对值的最大值 
			tmp = abs(Ab[c][r])
			if tmp > num_max:
				ind_max,num_max = c, tmp
		if num_max < epsilon:
			return None              # A为奇异矩阵，返回None
		Ab = swapRows(Ab,r,ind_max)  # 将绝对值最大值所在行交换到对角线元素所在行
		Ab = scaleRow(Ab,r,1/Ab[ind_max][r])  # 将列c的对角线元素缩放为1
		for c in range(col):         # 将列c的其他元素消为0
			if c == r:
				continue
			Ab = addScaledRow(Ab,c,r,-1*Ab[c][r]])
	Ab = matxRound(Ab, decPts)
	N = [0] * row
	for r in range(row):
		N[i] = Ab[r][col-1]
	
    return N
