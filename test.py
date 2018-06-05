import QuadraticSieveFactorization4 
import RhoFactorization
from random import randint
from QuadraticSieveFactorization4 import isprime

def makeprime(limit):
    while True:
        z = randint(limit - 100000, limit)
        if isprime(z):
            return z

#f = open('output.csv','w')
#print('bits\tnumber\tquadratic\trho', file = f)
#f.close()
#for j in range(30):
#    for i in range(10, 16):            
#        f = open('output.csv','a')
#        number = makeprime(10 ** i) * makeprime(10 ** i)
#        print('{}'.format(2 * i), file = f, end = '\t')
#        print('{}'.format(number), file = f, end = '\t')
#        print('{}'.format(QuadraticSieveFactorization4.test(number)), file = f, end = '\t')
#        print('{}'.format(RhoFactorization.test(number)), file = f) 
#        f.close()

for j in range(10):
    for i in range(16, 18):
        f = open('output.csv','a')
        number = makeprime(10 ** i) * makeprime(10 ** i)
        print('{}'.format(2 * i), file = f, end = '\t')
        print('{}'.format(number), file = f, end = '\t')
        print('{}'.format(QuadraticSieveFactorization4.test(number)), file = f, end = '\t')
        print('{}'.format(RhoFactorization.test(number)), file = f) 
        f.close()

for j in range(10):
    for i in range(18, 20):
        f = open('output.csv','a')
        number = makeprime(10 ** i) * makeprime(10 ** i)
        print('{}'.format(2 * i), file = f, end = '\t')
        print('{}'.format(number), file = f, end = '\t')
        print('{}'.format(QuadraticSieveFactorization4.test(number)), file = f, end = '\t')
        print('{}'.format(RhoFactorization.test(number)), file = f) 
        f.close()
