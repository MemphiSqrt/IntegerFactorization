pow_range = 4000
prime_range = 1000000


def gcd(x, y):
    return x if y == 0 else gcd(y, x % y)


def power(x, v):
    ans = 1
    while v != 0:
        if (v & 1) == 1:
            ans = ans * x
        x = x * x
        v >>= 1
    return ans


def prime_gan(n):
    prime = []
    p = [0] * n
    for i in range(2, n):
        if p[i] == 0:
            prime.append(i)
            p[i] = i
        for j in prime:
            if j > p[i] or i * j >= n:
                break
            p[i * j] = j
    return prime


prime_list = prime_gan(prime_range)


def prime_size(n):
    return prime_list[:n]


pow2 = [1] * pow_range
for i in range(1, pow_range):
    pow2[i] = pow2[i-1] * 2


def get_pow2(n):
    return pow2[n]
