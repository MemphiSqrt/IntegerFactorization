import math
from Util import prime_n_limit

n = int(input())
log_n = math.log(n)
B = math.ceil(math.e ** (0.56 * math.sqrt(log_n * math.log(log_n)))) + 300
print("B smooth is = {}".format(B))
print("prime_base is = {}".format(prime_n_limit(B)))
