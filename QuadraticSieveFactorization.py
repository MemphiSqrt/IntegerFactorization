from Util import *
from math import sqrt, ceil
from random import randint

# parameters
prime_check = prime_size(40)
prime_cnt = 1000
prime_elect_cnt = 400


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


def quadratic(num):
    factor = []

    # make prime and initialize
    prime = prime_size(prime_cnt)
    prime_set = [[]] * prime_cnt
    # initialize for quadratic residue
    for i in range(prime_cnt):
        prime_set[i] = [[] for _ in range(prime[i])]
        for j in range(prime[i]):
            prime_set[i][j ** 2 % prime[i]].append(j)

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
                if len(prime_set[i][n % prime[i]]) != 0:
                    # have quadratic residue
                    res = n % prime[i] # save mod value to accelerate program
                    if res == 0:
                        factor.append(prime[i])
                        quadratic_iteration(n // prime[i])
                        return
                    prime_table.append(prime[i])
                    tmp = set({positive_mod(prime_set[i][res][0] - sqrt_n, prime[i])})
                    # if res is not 0, there's two quadratic residues
                    if prime[i] != 2:
                        if len(prime_set[i][res]) < 2:
                            print('debug plz! error 4')
                            print('prime now is {}'.format(prime[i]))
                            exit(0)
                        tmp.update({positive_mod(prime_set[i][res][1] - sqrt_n, prime[i])})
                    quadratic_res.append(tmp)
                if len(prime_table) == prime_iter:
                    break

            if len(prime_table) < prime_iter:
                print('prime cnt is not enough!')
                print('prime iter = {}'.format(prime_iter))
                exit(0)

            print(prime_table)

            number_iter = 0
            while True:
                y_value = y_f(number_iter)
                if number_iter < 10:
                    print('number iter = {}, y_value = {}'.format(number_iter, y_value))

                pow2 = 1
                xor_vector = 0
                bit_vec = []
                for i in range(prime_iter):
                    if (number_iter % prime_table[i]) in quadratic_res[i]:
                        if y_value % prime_table[i] != 0:
                            print('debug plz! error 1 ')
                            exit(0)
                        else:
                            y_value //= prime_table[i]
                            xor_vector |= pow2
                        bit_vec.append(1)
                    else:
                        bit_vec.append(0)
                    pow2 <<= 1
                if y_value == 1:
                    print('prime_base cnt = {}, number iter = {}, count = {}'.format(prime_iter, number_iter, linear_inde[0]))
                    # print(sqrt_n + number_iter)
                    number_set = set({number_iter})
                    # if number_iter in [33, 326, 811, 267, 156, 511]: print('value = {}, iter = {}'.format(bit_vec, number_iter))
                    flag, result_set = linear_equation(xor_vector, number_set)
                    if flag:
                        if check(result_set):
                            return
                        else:
                            break

                number_iter += 1

    quadratic_iteration(num)
    return factor


number = randint(1, 10 ** 30)
# number = 74756495428169733753
print('number = {}'.format(number))
factor = quadratic(number)
print(factor)
print(number)
