from Util import *
from random import randint

# parameters
prime_check = prime_size(40)
prime_cnt = 50000
prime_elect_cnt = 5000
block_size = 1000000


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
    if n <= prime_check[len(prime_check) - 1]:
        return n in prime_check
    for i in prime_check:
        if witness(i, n):
            return False
    return True


def sqrt_int(n):
    k = 1
    prev1 = prev2 = k
    while abs(k * k - n) > 1:
        tmp = (k + n // k) // 2
        if tmp == prev1: break
        prev1 = prev2
        prev2 = k
        k = tmp
    return k


def quadratic_residue(a, p):
    if p == 2:
        return 1 if a == 1 else 0

    ch = power(a, (p - 1) // 2, p)
    if ch != 1: return -1
    for i in range(2, p):
        if power(i, (p - 1) // 2, p) == p - 1:
            b = i
            break
    p_1 = p - 1
    t = 0
    while p_1 % 2 == 0:
        t += 1
        p_1 //= 2
    s = p_1
    x = power(a, (s + 1) // 2, p)
    a_inv = power(a, p - 2, p)
    e = x * x % p * a_inv % p
    two_t = (p - 1) // s # 2^t
    two_i = 1
    for k in range(1, t):
        if power(e, two_t // two_i // 4, p) != 1:
            x = x * power(b, two_i * s, p) % p
        two_i *= 2
        e = x * x % p * a_inv % p
    return x


def quadratic(num):
    factor = []

    # make prime and initialize
    prime = prime_size(prime_cnt)

    def quadratic_iteration(n):
        if isprime(n):
            factor.append(n)
            return
        # print('work here')
        sqrt_n = sqrt_int(n)
        if sqrt_n * sqrt_n < n: sqrt_n += 1
        # print(sqrt_n * sqrt_n - n)
        print('n = {}, sqrt of n = {}'.format(n, sqrt_n))
        # exit(0)

        def y_f(x): return (x + sqrt_n) ** 2 - n

        def positive_mod(x, p): return (x % p + p) % p

        # equation_set = []
        # equation_vec = []
        def linear_equation(vec, index_set):
            pow2 = 1
            for i in range(prime_iter):
                if (vec & pow2) != 0:
                    if len(equation_set[i]) == 0:
                        linear_inde[0] += 1
                        equation_set[i] = index_set
                        equation_vec[i] = vec
                        return False, index_set
                    else:
                        vec ^= equation_vec[i]
                        index_set ^= equation_set[i]
                pow2 <<= 1
            if vec != 0:
                print('debug plz! error 2')
                exit(0)
            return True, index_set

        def check(index_set):
            x = 1
            y = 1
            # print(index_set)
            for i in index_set:
                x *= (sqrt_n + i)
                y *= y_f(i)
            sqrt_y = sqrt_int(y)
            # print('y = {}, sqrt_y = {}'.format(y, sqrt_y))
            # exit(0)
            if sqrt_y ** 2 != y:
                print('debug plz! error 3')
                print('y = {}, sqrt_y = {}'.format(y, sqrt_y))
                exit(0)
            y = sqrt_y
            if y > x:
                tmp = x
                x = y
                y = tmp
            z = gcd(x - y, n)
            if z != 1 and z != n:
                quadratic_iteration(n // z)
                quadratic_iteration(z)
                return True
            z = gcd(x + y, n)
            if z != 1 and z != n:
                quadratic_iteration(n // z)
                quadratic_iteration(z)
                return True
            return False

        for prime_iter in range(prime_elect_cnt, prime_elect_cnt + 1000):
            linear_inde = [0]
            equation_set = [set() for _ in range(prime_iter)]
            equation_vec = [[] for _ in range(prime_iter)]
            prime_table = []
            quadratic_res = []

            for i in range(prime_cnt):
                if n % prime[i] == 0:
                    factor.append(prime[i])
                    quadratic_iteration(n // prime[i])
                    return

            link = [[[] for _ in range(block_size)] for _ in range(2)]
            status = 0
            for i in range(prime_cnt):
                if prime[i] > block_size:
                    print('debug plz! error 5')
                    exit(0)

                res = n % prime[i]  # save mod value to accelerate program
                ans = quadratic_residue(res, prime[i])
                if ans != -1:
                    # print('work here')
                    it_ind = len(prime_table)
                    prime_table.append(prime[i])
                    uk = positive_mod(ans - sqrt_n, prime[i])
                    link[status][uk].append(it_ind)
                    if prime[i] != 2:
                        uk = positive_mod(prime[i] - ans - sqrt_n, prime[i])
                        link[status][uk].append(it_ind)
                if len(prime_table) == prime_iter:
                    break

            if len(prime_table) < prime_iter:
                print('prime cnt is not enough!')
                print('prime iter = {}'.format(prime_iter))
                exit(0)

            print(prime_table)

            number_iter = 0
            block_count = -1
            status = 1
            while True:
                if number_iter % block_size == 0:
                    # del link[status]    #  = [[] for _ in range(block_size)]
                    link[status] = [[] for _ in range(block_size)]
                    status ^= 1
                    block_count += 1
                y_value = y_f(number_iter)

                xor_vector = 0
                tmp_dis = block_count * block_size
                for i in link[status][number_iter - tmp_dis]:
                    p = prime_table[i]
                    if y_value % p != 0:
                        print('debug plz! error 1 ')
                        exit(0)
                    else:
                        y_value //= p
                        xor_vector |= get_pow2(i)
                    if number_iter + p >= tmp_dis + block_size:
                        link[status ^ 1][number_iter + p - tmp_dis - block_size].append(i)
                    else:
                        link[status][number_iter + p - tmp_dis].append(i)

                if y_value == 1:
                    print('prime_base cnt = {}, number iter = {}, count = {}'.format(prime_iter, number_iter, linear_inde[0]))
                    # print(sqrt_n + number_iter)
                    number_set = set({number_iter})
                    flag, result_set = linear_equation(xor_vector, number_set)
                    if flag:
                        if check(result_set):
                            return
                        else:
                            break

                number_iter += 1

    quadratic_iteration(num)
    return factor


number = randint(1, 10 ** 40)
# number = 15347
print('number = {}'.format(number))
factor = quadratic(number)
print(factor)
print(number)
