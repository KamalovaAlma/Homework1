import math
import sys

h=0
tau=0

B=0
A=0
C=0

n=5
m=40

alpha=[]
betta=[]

def main():
	global n
	global m
	for n in [5,15,25,35]:#
		#print('n=',n)

		for m in [40,80,160,320]:#
			h=math.pi/n
			tau=10/m

			u = [[0]*(m+1) for i in range(n+1)]

			u1 = execM1(u, h, tau)
			E = execE(u1, h, tau)
			print('M1: h=',h,' tau=',tau,' E=',E)
			#print(E)

			u = [[0]*(m+1) for i in range(n+1)]
			B=tau/(h**2)
			A=tau/(h**2)
			C=-(1 + (2*tau)/(h**2))
			u2 = execM2(u, A, B, C, h, tau)
			E = execE(u2, h, tau)
			print('M2: h=',h,' tau=',tau,' E=',E)
			#print(E)

			u = [[0]*(m+1) for i in range(n+1)]
			B_kn=tau/(2 * h**2)
			A_kn=tau/(2 * h**2)
			C_kn=-(1 + tau/(h**2))
			u3 = execM3(u, A_kn, B_kn, C_kn, h, tau)
			E = execE(u3, h, tau)
			print('M3: h=',h,' tau=',tau,' E=',E)
			#print(E)
		print('===========')
	

def execE(u, h, tau):
	E=0
	for i in range(0,n+1):
		for j in range(0,m+1):
			div = math.fabs(getUxt(h*i, tau*j) - u[i][j])
			if div>=E:
				E=div
	return E

def execM1(u, h, tau):
	for i in range(0, n + 1):
		u[i][0] = math.sin(i * h)

	for j in range(1, m + 1):
		tj=j*tau
		u[0][j] = math.log(tj**2 + 1)
		u[n][j] = math.log(tj**2 + 1)

	for j in range(0, m):
		for i in range(1, n):
			u[i][j+1]=(tau/h**2) * (u[i-1][j]-2*u[i][j]+u[i+1][j]) + tau*getGij(i*h, j*tau) + u[i][j]

	return u

def execM2(u, A, B, C, h, tau):
	#неявная схема
	alpha=[0]*(n + 1)
	betta=[0]*(n + 1)

	for i in range(0, n + 1):
		u[i][0] = math.sin(i * h)

	for j in range(1, m + 1):
		tj=j*tau
		u[0][j] = math.log(tj**2 + 1)
		u[n][j] = math.log(tj**2 + 1)	

	for j in range(0,m):
		alpha=[0]*(n + 1)
		betta=[0]*(n + 1)

		tj=j*tau

		alpha[1]=0
		betta[1]=math.log(tj ** 2 +1)

		for i in range(1,n):
			xi=i*h
			alpha[i+1] = -B/(A*alpha[i]+C)
			F=-u[i][j]-getGij(xi,tj)*tau
			betta[i+1]=(F-A*betta[i])/(A*alpha[i]+C)

		for i in reversed((range(1,n))):
			u[i][j+1] = alpha[i+1] * u[i+1][j+1] + betta[i+1]

	return u

def execM3(u, A, B, C, h, tau):
	#схема Кранка-Никольсона
	alpha=[0]*(n + 1)
	betta=[0]*(n + 1)

	for i in range(0, n + 1):
		u[i][0] = math.sin(i * h)

	for j in range(1, m + 1):
		tj=j*tau
		u[0][j] = math.log(tj**2 + 1)
		u[n][j] = math.log(tj**2 + 1)	

	for j in range(0,m):
		alpha=[0]*(n + 1)
		betta=[0]*(n + 1)

		tj=j*tau

		alpha[1]=0
		betta[1]=math.log(tj ** 2 +1)

		for i in range(1,n):
			alpha[i+1] = -B/(A*alpha[i]+C)
			F=-(tau*u[i-1][j]/(2*h**2) + (1-tau/(h**2))*u[i][j] + tau*u[i+1][j]/(2*h**2) + tau*getGij(i*h, j*tau))
			betta[i+1]=(F-A*betta[i])/(A*alpha[i]+C)

		for i in reversed((range(1,n))):
			u[i][j+1] = alpha[i+1] * u[i+1][j+1] + betta[i+1]

	return u


def printU(u):
	for i in range(0,n+1):
		print(u[i])	

#u(x; t) = sin(x) + ln(t^2 + 1):
def getUxt(x,t):
	return math.sin(x) + math.log(t ** 2 +1)

def getGij(xi, tj):
	return math.sin(xi) + (2*tj)/(tj**2 + 1)


if __name__ == "__main__": main()