from gj_solution import gj_Solve

def calculateMSE(X,Y,m,b):
    ind = len(X)
    if ind == 0:
        return "point sets error"
    mse = [(y-m*x -b)**2 for x,y in zip(X,Y)]
    calmse = sum(mse)/ind
    return calmse


#3.4 XtXh = XtY

def linearRegression(X, Y):
    leng = len(X)
    XtX = [[0,0],[0,0]]
    XtY = [[0],[0]]
    XtX[0][0] = sum([x**2 for x in X])
    XtX[0][1] = XtX[1][0] = sum([x for x in X])
    XtX[1][1] = leng
    XtY = [[sum([x*y for x,y in zip(X,Y)])],
           [sum([y for y in Y])]]
    m,b = gj_Solve(XtX, XtY)
    return m,b

m2,b2 = linearRegression(X,Y)
assert isinstance(m2,float),"m is not a float"
assert isinstance(b2,float),"b is not a float"
print(m2,b2)
