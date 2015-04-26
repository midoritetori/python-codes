# -*- coding: utf-8 -*-

"""
@author: tetori
tetori@naver.com

Short instruction:

Structure of p-adic numbers works like as the floating point.
It divides two parts: significand and exponent. (Base = prime.)
Significand is an element of (Z/prime**expmax*Z).
(Remember the algebraic definition of p-adic integer.)

expmax = the maximal length of significand.
padic.digit = significand of padic.
padic.pexp = exponent of padic.

How to input a p-adic number:
if you enter
>> x = padic(digit, pexp)
then x == digit * prime**pexp. If you does not give pexp,
then pexp is automatically setted to 0.

abs(x) : calculate p-adic absolute value of x
ordinary arithematical operation such as addition, multiplication
works well.
x.inv() calculates the multiplicative inverse of x.
x**y is defined only when y is int or long.

sqrt(x): calculate p-adic square root. It works when the
square root of x exists.

exp(x): calculcate p-adic exponentational. 
log(x): calculate p-adic logarithm.

"""
prime = 5 # Setting the prime you want!
expmax = 30 # Maximal digit calculated.
reprtype = "array"

"""
reprtype list:
series: "a(-k) * p^(-k) + ... + a(0) * p^0 + a(1) * p^1 + ..."
(each a(i) is an integer lies between 0 and p-1)

array: "a(-k), ... , "." , a(0), a(1), ..."

float: "number * p**pexp"

primitive: (digit, pexp)

"""
import math

ub = prime**int(expmax*1.125) # upper bound of the digit.
trueub = prime**expmax

class padic(object):
    # Class padic is immutable. 
    def __setattr__(sel, *args):
        raise TypeError("Can not modify immutable instance")
    __delattr__ = __setattr__
    
    def __init__(self, digit, pexp = 0):
        if(digit == 0):
            pexp = 0
        elif(digit < 0):
            digit = ub + digit
        
        while(not digit % prime and digit):
            digit, pexp = digit // prime, pexp+1
        while(digit >= ub):
            digit %= ub
        
        super(padic, self).__setattr__('digit', digit)
        super(padic, self).__setattr__('pexp', pexp)
        
    def __add__(a, b):
        if type(b) in [int, long]:
            b = padic(b)
        if(a.pexp>=b.pexp):
            resultdigit = (a.digit * prime**(a.pexp - b.pexp) + b.digit) % ub
            resultpexp = b.pexp
            return padic(resultdigit, resultpexp)
        else:
            resultdigit = (b.digit * prime**(b.pexp - a.pexp) + a.digit) % ub
            resultpexp = a.pexp
            return padic(resultdigit, resultpexp)
            
    def __mul__(a, b):
        if type(b) in [int, long]:
            b = padic(b)
        resultdigit = (a.digit*b.digit) % ub
        return padic(resultdigit, a.pexp+b.pexp)
        
    def __sub__(a, b):
        if type(b) in [int, long]:
            b = padic(b)
        if(a.pexp>=b.pexp):
            resultdigit = (a.digit * prime**(a.pexp - b.pexp) - b.digit) % ub
            resultpexp = b.pexp
            return padic(resultdigit, resultpexp)
        else:
            resultdigit = (-b.digit * prime**(b.pexp - a.pexp) + a.digit) % ub
            resultpexp = a.pexp
            return padic(resultdigit, resultpexp)  
    
    def inv(self):
        if(self.digit == 0):
            raise(ZeroDivisionError("%d-adic division by zero" % prime))
        
        # It uses extended Euclid algorithm
        r0, r1 = trueub, self.digit
        t0, t1 = 0, 1
        while(r1):
            if(r0 < r1):
                r0 , r1 , t0 , t1 = r1, r0, t1, t0
            q1 = r0 // r1
            r0, r1 = r1, r0 % r1
            t0, t1 = t1, t0 - q1 * t1
        
        return padic(t0 % ub, -self.pexp)
    
    def __div__(a, b):
        if type(a) in [int, long]:
            a = padic(a)
        if type(b) in [int, long]:
            b = padic(b)
        return padic.__mul__(a, b.inv())
    
    def __pow__(a, b):
        if(type(b) != int and type(b) != long):
            raise(TypeError("Exponent of p-adic exponentation must be an integer."))
        
        if(b < 0):
            a, b = padic.__rdiv__(a,1), -b
        elif(b == 0):
            return padic(1)
        elif(a.digit == 1):
            return padic(1, a.pexp*b)
        
        bsave = b
        rev_bsave = 1
        result = 1
        while(bsave):
            rev_bsave <<= 2
            rev_bsave += 3 if (bsave%2) else 2
            bsave >>= 1
        while(rev_bsave > 1):
            result = (result**2)%ub if (rev_bsave%4 == 2) else (result**2 * a.digit)%ub
            rev_bsave >>= 2
        
        return padic(result % ub, a.pexp*b)
        
    def __neg__(a):
        if(a.digit != 0):
            return padic(ub-a.digit, a.pexp)
        else:
            return padic(0,0)
    
    def __abs__(a):
        return padic(1,-a.pexp)
        
    def __eq__(a, b):
        if(type(a) in [int, long]):
            return (a == b.digit * (prime ** b.pexp) )
        elif(type(b) in [int, long]):
            return (b == a.digit * (prime ** a.pexp) )
        else:
            return (a.digit == b.digit) and (a.pexp == b.pexp)
    
    def __ne__(a, b):
        if(type(a) in [int, long]):
            return (a != b.digit * (prime ** b.pexp) )
        elif(type(b) in [int, long]):
            return (b != a.digit * (prime ** a.pexp) )
        else:
            return (a.digit != b.digit) or (a.pexp != b.pexp)
    
    def __pos__(a):
        return a
        
    def __radd__(a,b):
        if type(b) in [int, long]:
            b = padic(b)
        return padic.__add__(b,a)
        
    def __rmul__(a,b):
        if type(b) in [int, long]:
            b = padic(b)
        return padic.__mul__(b,a)
        
    def __rsub__(a,b):
        if type(b) in [int, long]:
            b = padic(b)
        return padic.__sub__(b,a)
        
    def __rdiv__(a,b):
        if type(b) in [int, long]:
            b = padic(b)
        return padic.__div__(b,a)
        
    def __rpow__(a,b):
        return padic.__pow__(b,a)
    
    def __int__(a):
        return int(a.digit * prime**a.pexp)
    
    def __float__(a):
        return a.digit * prime**a.pexp
    
    """
    intype : representing stype. default is reprtype.
    inexpmax : maximum representation digit. default is expmax. 
               this variable only needed for "series" and "array"
    omitzero : (in "series" type) represent zero-coefficient terms lies between
               two nonzero terms. default is False.
    """
    def __repr__(self, intype = reprtype, inexpmax = expmax, omitzero = False):
        if(intype == "series"):
            digit = self.digit
            i = 0
            series = ""
            if(digit == 0):
                return "0"
            
            while(digit and i < inexpmax ):
                if(digit%prime != 0 and omitzero):
                    series = series + ("%d * %d^%d + " % ((digit % prime, prime, i+self.pexp)))
                digit, i = digit//prime, i+1
            return series[:-3]
            
        elif(intype == "array"):
            digit = self.digit
            i = 0
            array = [0] * self.pexp
            if(digit == 0):
                return "0"
            
            while(digit and i < inexpmax):
                array.append(int(digit % prime))
                if(i+self.pexp == -1):
                    array.append(".")
                digit, i = digit//prime, i+1
            return str(array)[1:-1]
            
        elif(intype == "number"):
            return "%d * %d^%d" % (self.digit, prime, self.pexp)
            
        elif(intype == "primitive"):
            return "(%d,%d)" % (self.digit, self.pexp)
        
        else:
            raise(TypeError("variable 'reprtype' is invalid."))
    
    def __str__(self, intype = reprtype, inexpmax = expmax):
        return padic.__repr__(self, intype, inexpmax = expmax)
        
def inv(x):
    if type(x) in [int, long]:
        x = padic(x)
    return x.inv()
    
def JacobiSymb(a, b):
    if(b%2 == 0):
        raise(ValueError("Second argument of Jacobi symbol should be odd."))
    
    value = 1
    while(a >= 2):
        if(a%2 == 0):
            a, value = a/2, value * (1 if (b%8 == 1 or b%8 ==7) else (-1))
        if(a>b):
            a %= b
        if(a<b and a%2 == 1):
            a, b = b, a
            if(a%4 == 3 and b%4 == 3):
                value *= -1
        
    return value

# Calculate the least quadratic residue of x.num % prime modulo prime
# by using Tonelli–Shanks algorithm.
def quadres(n):
    if(JacobiSymb(n, prime) == -1):
        raise(ValueError("Quadratic residue of given argument does not exist"))
    
    varQ = (prime - 1)/2
    varS = 1
    while(not varQ % 2):
        varQ, varS = varQ/2, varS+1

    if(varS == 1):
        varR = (n**((prime+1)/4)) % prime
    
    else:
        nonres = 2
        while(JacobiSymb(nonres, prime) == 1):
            nonres += 1
        
        varC = (nonres ** varQ) % prime
        varR, varT, varM = (n ** ((varQ+1)/2) )%prime, (n**varQ)%prime, varS
        
        while(varT != 1):
            i = 0
            temp = varT
            while(i < varM and temp != 1):
                temp = (temp**2)%prime
                i += 1
            
            varB = (varC**(2**(varM - i - 1) ) )%prime
            varR = (varR * varB)%prime
            varT = (varT * varB ** 2)%prime
            varC, varM = (varB**2)%prime, i
        
    if(varR > prime/2):
        varR = prime - varR
        
    return varR

def sqrt(x):
    if type(x) in [int, long]:
        x = padic(x)
    
    if(x.digit == 0):
        return 0
        
    if(x.pexp % 2):
        raise(ValueError("sqrt of the given value does not exist"))

    if(prime == 2):
        seed = x.digit % 8
        if(seed == 1):
            a = 1
        else:
            raise(ValueError("sqrt of the given value does not exist"))
    
    else:
        seed = x.digit % prime
        try:
            a = quadres(seed)
        except(ValueError):
            raise(ValueError("sqrt of the given value does not exist"))
    
    j = 0 if (prime-2) else 2
    while(j <= expmax):
        try:
            a, j = a - (a*a - x.digit)*inv(a+a), j+1
        except(NameError):
            raise(ValueError("sqrt of the given value does not exist"))
    return padic(a.digit, x.pexp/2)
    
def exp(x):
    if type(x) in [int, long]:
        x = padic(x)
    if(x.pexp < (1 if (prime>2) else 2)):
        raise(ValueError("exponentation of given value cannot be computable"))
        
    """
    Set a summation range. Let nu(p,x) be a p-valuation of 
    p-adic number x. Then we get
    
        nu(p, x^n/factorial(n)) >= n*nu(p,x) + (S(n)-n)/(p-1)
    
    where S(n) is the sum of digits of n written in base p.
    (See Katok's lecture note about p-adic analysis written in 2001.)
    We shall make nu(p, x^n/factorial(n)) > expmax + nu(p,x).
    (The reason why we add nu(p,x) to expmax is we are using
     float-point data type.)
    From S(n) >= 1 for all n >= 1, we get
    
        n >= [(p-1)*(expmax + nu(p,x)) -1] / [(p-1)*nu(p,x) - 1]
    
    """
    
    sumbound = ((prime-1) * (expmax+x.pexp) - 1) / float((prime-1)*x.pexp - 1)
    value = 1
    term = 1
    k = 1
    while(k <= sumbound):
        term *= x/k
        value += term
        k += 1
    
    return value


def LambertWm1(x0):
    if(x0 >= 0 or x0 < -1/math.e):
        raise(ValueError("Lambert W function is not defined for given value"))
    w = -1
    error = 1
    n = 1
    """
    We will use following iteration:
    w(0) = -1, w(n+1) = log(-x0)- log(-w(n))
    
    To estimate the error, we will derive following inequality:
    |w(n+1) - w(n)|  = |log(-w(n)) - log(-w(n-1))|
                    <= |log(1 + (w(n) - w(n-1) )/w(n-1) )|
                    <= |w(n) - w(n-1)| / |w(n-1)|
    
    From this we get
    |w(n+1) - w(n)| <= |w(1) - w(0)| / |w(0) * w(1)^(n-1)|
    
    Therefore
    |w(inf) - w(n)| <= sum k from n to inf |w(k+1) - w(k)|
                    <= |w(1) - w(0)| / |w(0) * w(1)^(n-2) * (|w(1)| - 1)|
    """
    
    while(abs(error) > 10**(-15)):
        w = math.log(-x0) - math.log(-w)
        n += 1
        error = abs( (1 - math.log(-x0)) / (math.log(-x0) ** (n-2) * (1 + math.log(-x0))) )
        
    return w


def log(x):
    if type(x) in [int, long]:
        x = padic(x)
    x1 = x-1
    if(x1==0):
        return 0
    if(x1.pexp < 1):
        raise(ValueError("logarithm of given value is not computable or may not exist"))
    
    """
    Set a summation range. If n is a max. bound then
    
        nu(p, x^n / n) = n*nu(p,x) - nu(p,n)
        
    We shall make nu(p, x^n / n) >= expmax.
    From the inequality nu(p, x) <= log(x)/log(p),
    we consider the inequality
    
        n*nu(p,x) - log(n)/log(p) >= expmax + nu(p,x).
    
    Long and tedious calculation gives a bound:
    
        n >= -1/(nu(p,x)*log(p)) * W(-nu(p,x)*log(p)*p^(-expmax-nu(p,x)) )
    
    where W is a Lambert W-function defined over (-1/e,0) (i.e. W_{-1}.)
    
    """
    sumbound = 1 - 1/(x1.pexp * math.log(prime)) * LambertWm1(-x1.pexp * math.log(prime) * prime**(- expmax - x.pexp))
    # sumbound algorithm shall be given
    value = padic(0)
    k = 1
    while(k < sumbound):
        value += (-1)**(k-1) * x1**k/k
        k += 1
    
    return value
    
# not supported
# def polysol(*pargs):
#    pass
