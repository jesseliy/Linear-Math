# TODO 实现 Gaussain Jordan 方法求解 Ax = b

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

def swapRows(M, r1, r2):
    M[r1],M[r2] = M[r2],M[r1]
    pass

def scaleRow(M, r, scale):
    if scale == 0:
        raise ValueError
    else:
        M[r] = [x*scale for x in M[r]]
    pass

def addScaledRow(M, r1, r2, scale):
    M[r1] = [x+y*scale for x,y in zip(M[r1],M[r2])]
    pass

def augmentMatrix(A, b):
    N = [x+y for x,y in zip(A,b)]
    return N

def gj_Solve(A, b, decPts=4, epsilon = 1.0e-16):

    row = len(A[0])

    if row != len(b):
        return None

    Ab = augmentMatrix(A,b)

    for i in range(len(Ab)):



    return None
