def gcd(x, y):
    return x if y == 0 else gcd(y, x % y)


def prime_gan(n):
    prime = []
    p = [[0] * n]
    for i in range(2, n):
        if p[i] == 0:
            prime.append(i)
            p[i] = i
        for j in prime:
            if j > p[i]:
                break
            p[i * j] = j
    return prime


prime_list = prime_gan(1000000)


def prime_size(n):
    return prime_list[:n]
