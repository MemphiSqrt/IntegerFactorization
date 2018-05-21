from QuadraticSieveFactorization4 import quadratic
from QuadraticSieveFactorization4 import factor
from time import time

# number = randint(1, 10 ** 30)
number = 777796089990233610621609434355
print('number = {}'.format(number))
time_cnt = time()
quadratic(number)
print(factor)
print('time = {}'.format(time() - time_cnt))
print(number)