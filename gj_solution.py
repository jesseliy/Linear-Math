
def transpose(M):
    N = []
    row = len(M)
    col = len(M[0])
    for i in range(col):
        tmp = []
        for j in range(row):
            tmp.append(M[j][i])
        N.append(tmp)
    return N


def matxMultiply(A, B):
    row = len(A)
    col = len(A[0])
    N = []
    if not(len(B)== col):
        raise ValueError
    else:
        BT = transpose(B)
        for i in range(row):
            tmp = []
            for j in range(len(B[0])):
                V_tmp = [x*y for x,y in zip(A[i], BT[j])]
                tmp.append(sum(V_tmp))
            N.append(tmp)
    return N


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
    if len(A[0]) != len(b):
        return None
    Ab = augmentMatrix(A,b)
    row = len(Ab)
    col = len(Ab[0])
    for c in range(col-1):
        ind_max = c
        num_max = abs(Ab[c][c])
        for r in range(c+1,row):
            tmp = abs(Ab[r][c])
            if tmp > num_max:
                ind_max,num_max = r, tmp
        if num_max < epsilon:
            return None
        Ab = swapRows(Ab,c,ind_max)
        lam = 1.0/Ab[c][c]
        Ab = scaleRow(Ab,c,lam)
        for r in range(row):
            if r == c:
                continue
            Ab = addScaledRow(Ab,r,c,-1*Ab[r][c])
    Ab = matxRound(Ab, decPts)
    N = [0] * row
    for r in range(row):
        N[r] = Ab[r][col-1]
    return N
