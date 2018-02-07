import backgroundcode
import gj_Solve from gj_solution

# 3.2.2 calculate MSE of a line
'''
input X,Y: Points
input m,b: Line's parameter: y = mx + b (float)
'''  
 
def calculateMSE(X,Y,m,b):
    ind = len(X)
	if ind == 0:
    	return "point sets error"
    K = [(y-m*x -b)**2 for x,y in zip(X,Y)]
    return sum(K)/ind

print(calculateMSE(X,Y,m1,b1))



#3.4 XtXh = XtY

'''
input X,Y: Points
input m,b: Line's parameter: y = mx + b
h = [m,b]
'''
def linearRegression(X, Y):
    leng = len(X)
	XtX = [0,0] * 2
	XtY = [0] * 2
	XtX[0][0] = sum([x**2 for x in X])
	XtX[0][1] = XtX[1][0] = sum([x for x in X])
	XtX[1][1] = leng
	XtY = [sum([x*y for x,y in zip(X,Y)]),
	       sum([y for y in Y])]
	m,b = gj_Solve(XtX, XtY)

    return m,b

m2,b2 = linearRegression(X,Y)
assert isinstance(m2,float),"m is not a float"
assert isinstance(b2,float),"b is not a float"
print(m2,b2)