# background codes

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


	
