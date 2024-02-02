import math,random,time
from cmath import *
import decimal as dc


from itertools import groupby

g_gamma = 7
p_gamma = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]

def gamma(z):
	z = complex(z)

	if z.real < 0.5:
		return pi / (sin(pi*z)*gamma(1-z))
	else:
		z -= 1
		x = p_gamma[0]
		for i in range(1, g_gamma+2):
			x += p_gamma[i]/(z+i)
		t = z + g_gamma + 0.5

		return sqrt(2*pi) * t**(z.real+0.5) * exp(-t.real) * x

def incomplete_gamma(s,x,alternate=False):
	#Described in computing the incomplete Gamma function to arbitrary precision  - Serge Winitzki
	summer = 0

	if alternate:
		fast = False
		if s > 0:
			for n in range(0,150):
				try:
					if fast:
						#using log - not as accurate
						num = log((x**(n+s))).real
						den = factorialRamanujan(n,logged=True)
						bit = (1.0/(s+n))*((-1)**n)*(e**(num - den))
					else:
						num = ((-1)**n)*(x**(n+s)).real
						den = ((s + n)*factorial(n))
						bit = num/den

					try:
						summer += bit
					except:
						print(("Error", s,x,num,den))
				except:
					pass

		return (gamma(s) - summer).real
	else:
		i = factorial(s-1)
		j = e**(-x)
		l = 0
		for k in range(0,s):
			l +=  float((x**k).real)/factorial(k)

		return  (i*j*l).real


def factorialRamanujan(m,s=0,logged=True):
	if m == 0:
		if logged:
			return log(1)
		else:
			return 1


	else:
		if logged:
			return (m*log(m)) - m + (log(m*(1+4*m*(1+2*m)))/6) + (log(pi)/2).real
		else:
			return int(math.ceil(e**((m*log(m)) - m + (log(m*(1+4*m*(1+2*m)))/6) + (log(pi)/2)).real))

def factorial(m,s=0):
	value = int(1)
	if m != 0:
		while m != s:
			value = value*m
			m = m - 1
	return value


def uniform_product_density(p,n):
	a = (-1)**(n-1)
	b = factorialRamanujan(n-1)
	c = log(p).real**(n-1)

	return (a/b)*c

def cum_uniform_product(p,n,alternate = False):
	try:
		if p == 1:
			return 1
		else:
			if alternate:

				a = (-1)**(n)
				b = incomplete_gamma(n,-log(p))
				c = (log(p).real)**n
				d = factorialRamanujan(n-1,logged=False)
				e = (-log(p).real)**n
				return (a*b*c)/(d*e)

			else:
				a = incomplete_gamma(n,-log(p))
				b = factorial(n-1)
				return a/b
			
	except Exception as e:
		print(("Error in Sig correction:",p,n))
		print(e)
		return -1

def erf(z):
	summer = 0
	for i in range(0,100):

		try:
			num = ((-1.0)**i) * z**(2*i+1)
		except OverflowError:
			print("Overflow")

		denum = factorial(i)*(2*i+1)
		summer += num/denum


	erf = summer*(2 /math.sqrt(math.pi))
	if erf > 1 or erf < -1:
		erf = 1

	return erf


def product(values):
	prod = dc.Decimal(str(1))
	for v in values:
		prod *= dc.Decimal(str(v))

	return prod


def gini_coefficient(column):
	diffsum = 0
	try:
		for i, xi in enumerate(column[:-1], 1):
			diffsum += sum([abs(xi -y) for y in column[i:]])

		return  diffsum / (len(column)**2 * sum(column)/len(column))
	except:
		return 0

def mean_absolute_error(column_a, column_b):
	mae = sum([abs(column_a[i] - column_b[i]) for i in range(0,len(column_a))])/(sum(column_a)+sum(column_b))
	return mae

def pearson_correlation(column_a, column_b):
	try:
		products_mean = sum([column_a[i] * column_b[i] for i in range(0, len(column_a))])/len(column_a)
		covariance = products_mean - (sum(column_a)/len(column_a) * sum(column_b)/len(column_b))
		column_a_standard_deviation =  standard_deviation(column_a)
		column_b_standard_deviation = standard_deviation(column_b)
		pearson_correlation = covariance / (column_a_standard_deviation * column_b_standard_deviation)
		return pearson_correlation
	except:
		return 0

def standard_deviation(column):
	column_squared = [i**2 for i in column]
	squares_mean = sum(column_squared)/len(column_squared)
	column_mean = sum(column)/len(column)
	column_mean_squared = column_mean**2
	variance = squares_mean - column_mean_squared
	std_dev = math.sqrt(variance)
	return std_dev 

##---------------------------------##