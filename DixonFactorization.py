from Util import *
from math import sqrt, ceil
import numpy as np
from random import randint

def dixon(n):
    for prime_base in range(4, 1000):
        prime_table = prime_size(prime_base)
        xor_table = [-1] * prime_base
        xor_labels = [set([])] * prime_base
        store_vec_table = []
        store_xi_table = []

        def number2vector(x):
            vec = []
            xor_vec = 0
            if x == 0: return vec, -1, False;
            for i in prime_table:
                counter = 0
                while x % i == 0:
                    counter += 1
                    x //= i
                vec.append(counter)
                xor_vec = (xor_vec << 1) | (counter & 1)
            return vec, xor_vec, x == 1
            # return vector, xor vector, if prime_base enough

        def vector2number(vec):
            ans = 1
            for i in range(prime_base):
                ans = ans * power(prime_table[i], vec[i])
            return ans

        def check_list(xi, se):
            vec, xor_vec, flag = se
            if not flag:
                return False, set([])
            store_vec_table.append(vec)
            store_xi_table.append(xi)
            labels = set([len(store_xi_table) - 1])
            for i in range(prime_base):
                if (xor_vec & get_pow2(prime_base - i - 1)) != 0:
                    if xor_table[i] == -1:
                        xor_table[i] = xor_vec
                        xor_labels[i] = labels
                        break
                    else:
                        xor_vec ^= xor_table[i]
                        labels = (labels - xor_labels[i]) | (xor_labels[i] - labels)
            if xor_vec == 0:
                return True, labels
            return False, set([])

        def check_labels(labels):
            vec_sum = [0] * prime_base
            a = 1
            for ele in labels:
                vec_sum = list(np.array(vec_sum) + np.array(store_vec_table[ele]))
                a *= store_xi_table[ele]
            w = vector2number(vec_sum)

            z = gcd(a - w, n)
            if z != n and z != 1:
                return z
            z = gcd(a + w, n)
            if z != n and z != 1:
                return z
            return 1

        # in_number = ceil(sqrt(n))
        while True:
            # print(in_number)
            in_number = randint(1, n)
            print(in_number)
            flag, labels = check_list(in_number * in_number, number2vector(in_number * in_number % n))
            if flag:
                z = check_labels(labels)
                if z != 1:
                    return z
                break
            # in_number += 1


factor = dixon(1000000016000000063)
print(factor)
