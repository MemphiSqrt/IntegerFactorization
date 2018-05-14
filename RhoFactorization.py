from Util import prime_size
from random import randint
import time

prime = prime_size(15)


def gcd(a, b):
    return a if b == 0 else gcd(b, a % b)


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


print(prime)
start_time = time.time()
number = randint(1, 10 ** 40)
print(number)
print(rho(number))
end_time = time.time()
print(end_time - start_time)
