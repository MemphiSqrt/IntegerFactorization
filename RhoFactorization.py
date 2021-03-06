from Util import *
from random import randint
import time

prime = prime_size(20)


def power(x, v, p):
    ans = 1
    while v != 0:
        if (v & 1) != 0:
            ans = ans * x % p
        x = x * x % p
        v >>= 1
    return ans


def witness(a, n):
    u = n - 1
    t = 0
    while u % 2 == 0:
        t += 1
        u >>= 1
    x = power(a, u, n)
    for i in range(t):
        las = x
        x = x * x % n
        if x == 1 and las != 1 and las != n - 1:
            return True
    return x != 1


def isprime(n):
    if n <= prime[len(prime) - 1]:
        return n in prime
    for i in prime:
        if witness(i, n):
            return False
    return True


def rho(num):
    factor = []

    def rho_iteration(n):
        if isprime(n):
            factor.append(n)
            return
        x = xk = randint(1, n)
        cnt = cntdown = 2
        c = randint(1, n + 1)
        while True:
            cntdown -= 1
            x = x * x % n + c
            x = x if x < n else x - n
            if xk == x:
                xk = x = randint(1, n)
                cnt = cntdown = 2
                c = randint(1, n + 1)
                continue
            t = x - xk + n
            t = t if t < n else t - n
            t = gcd(t, n)
            if t > 1 and t < n:
                rho_iteration(t)
                rho_iteration(n // t)
                return
            if cntdown == 0:
                cnt <<= 1
                cntdown = cnt
                xk = x

    rho_iteration(num)
    return factor
    

def makeprime(limit):
    while True:
        z = randint(limit - 100000, limit)
        if isprime(z):
            return z


def test(number):
    print('number = {}'.format(number))
    start_time = time.time()
    print(rho(number))
    end_time = time.time()
    print('time = {}'.format(end_time - start_time))
    timec = end_time - start_time
    return timec
    
