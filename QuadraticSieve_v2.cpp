#include <iostream>
#include <gmpxx.h>
#include <bitset>
#include <vector>
#include <gmp.h>
#include <ctime>
#include <cmath>
#include <thread>
using namespace std;
typedef mpz_class Int;
#define rep(i,x,y) for(int i = x;i <= y; ++i)
#define dep(i,x,y) for(int i = x;i >= y; --i)

const int PRIME_CHECK_CNT = 40;
const int PRIME_LIM = 3000;
const int PRIME_CNT = 60000;
const int PRIME_UPP = 1000000;
const int SIEVE_CHUNK = 2000000;

double temper[SIEVE_CHUNK];
int prime[PRIME_CNT], p[PRIME_UPP], nip[SIEVE_CHUNK];
vector<int> link[SIEVE_CHUNK];
//bool ele_prime[SIEVE_CHUNK];

Int power(Int x, Int v, Int p) {
	Int ans = 1;
	for(; v != 0; v >>= 1, x = x * x % p)
		if((v & 1) == 1)
			ans = ans * x % p;
	return ans;
}

int power(int x, int v, int p) {
	int ans = 1;
	for(; v != 0; v >>= 1, x = (long long) x * x % p)
		if((v & 1) == 1)
			ans = (long long) ans * x % p;
	return ans;
}

Int gcd(Int a, Int b) {
	return b==0 ? a : gcd(b, a % b);
}

void prime_gen() {
	int cnt = 0;
	rep(i, 2, PRIME_UPP - 1) {
		if (p[i] == 0) {
			p[i] = i;
			prime[cnt++] = i;
			if (cnt == PRIME_CNT) return;
		}
		for (int j = 0; j < cnt && i * prime[j] < PRIME_UPP && p[i] >= prime[j]; ++j) {
			p[i * prime[j]] = prime[j];
		}
	}
	printf("PRIME_CNT if not enough 1!");
	exit(0);
}

bool witness(int a, Int n) {
	Int u = n - 1;
	int t = 0;
	while ((u & 1) == 0) {
		t++;
		u >>= 1;
	}
	Int x = power(a, u, n);
	Int las;
	rep(i, 1, t) {
		las = x;
		x = x * x % n;
		if (x == 1 && las != 1 && las != n - 1)
			return true;
	}
	return x != 1;
}

bool isprime(Int n) {
	rep(i, 0, PRIME_CHECK_CNT - 1)
    	if (n == prime[i])
			return true;
	if (n <= prime[PRIME_CHECK_CNT - 1])
		return false;
	rep(i, 0, PRIME_CHECK_CNT - 1)
		if (witness(prime[i], n))
			return false;
	return true;
}

int quadratic_residue(int a, int p) {
    if (p == 2)
        return a==1 ? 1 : 0;

    int ch = power(a, (p - 1) / 2, p);
    if (ch != 1) return -1;
    int b = 0;
    rep(i, 2, p-1)
        if (power(i, (p - 1) / 2, p) == p - 1) {
			b = i;
            break;
        }
    int p_1 = p - 1;
    int t = 0;
    while (!(p_1 & 1)) {
        t++;
        p_1 /= 2;
    }
    int s = p_1;
    int x = power(a, (s + 1) / 2, p);
    int a_inv = power(a, p - 2, p);
    int e = (long long) x * x % p * a_inv % p;
    int two_t = (p - 1) / s;
    int two_i = 1;
    rep(k, 1, t - 1) {
        if (power(e, two_t / two_i / 4, p) != 1)
            x = (long long) x * power(b, two_i * s, p) % p;
        two_i *= 2;
        e = (long long) x * x % p * a_inv % p;
    }
    return x;
}

inline Int y_f(const Int &x, const Int &sqrt_n, const Int &n) {
	Int Z = x + sqrt_n;
	return Z * Z - n;
}

unsigned int toInt(const Int &x) {
	unsigned int b = 1, ans = 0;
	rep(i, 1, 32){
		if ((b & x) != 0)
		ans |= b;
		b <<= 1;
	}
	//if(x != ans) puts("???");
	return ans;
}

double log2(Int n) {
	Int z = 1;
	int ans = 0;
	while (n>z) {
		z <<= 1;
		ans += 1;
	}
	return ans;
}

int positive_mod(const Int &x, int p) {
    return toInt((x % p + p) % p);
}

vector<Int> ans;
bitset<PRIME_LIM> A[PRIME_LIM], idA[PRIME_LIM], Ax, idAx, T;
int id[PRIME_LIM];
int A_cnt = 0;

void quadratic_sieve(Int n);

bool check(const Int &sqrt_n, const Int &n, const bitset<PRIME_LIM> &C) {
	Int X = 1;
	Int Y = 1;
	rep(i, 0, PRIME_LIM - 1) if (C[i]) {
		Int z = (sqrt_n + id[i]);
	//	cout<<z<<endl;
		X *= z;
		z = z * z;
		Y *= z - n;
	}
	Y = sqrt(Y);
	//cout<<X<<" "<<Y<<endl;
	//exit(0);
	Int g = gcd(X - Y, n);
	if (g != 1 && g != n) {
		quadratic_sieve(g);
		quadratic_sieve(n / g);
		return true;
	}
	g = gcd(X + Y, n);
	if (g != 1 && g != n) {
		quadratic_sieve(g);
		quadratic_sieve(n / g);
		return true;
	}
	return false;
}

int TT_clock = 0;

void quadratic_sieve(Int n) {
	if (isprime(n)) {
		ans.push_back(n);
		return;
	}
	rep(i, 0, PRIME_CNT - 1) {
		if (n % prime[i] == 0) {
			ans.push_back(prime[i]);
			quadratic_sieve(n / prime[i]);
			return;
		}
	}

	int list_size = 0;
	vector<int> list_prime, list_res;
	bool flag_prime = false;
	rep(i, 0, PRIME_CNT - 1) {
		int modn = toInt(n % prime[i]);
		int res = quadratic_residue(modn, prime[i]);
		if (res != -1) {
			list_prime.push_back(prime[i]);
			list_res.push_back(res);
			list_size++;
		}
		if (list_size == PRIME_LIM) {
			flag_prime = true;
			break;
		}
	}
	if(!flag_prime) {
		puts("PRIME_CNT is not enough!");
		exit(0);
	}

	rep(i, 0, PRIME_LIM - 1) A[i].reset(), idA[i].reset();
	A_cnt = 0;
	Int sqrt_n = sqrt(n);
	if(sqrt_n * sqrt_n < n) sqrt_n++;

	//rep(i,0,list_size-1) printf("%d ", list_prime[i]); puts("");
	//rep(i,0,list_size-1) printf("%d ", positive_mod(list_res[i] - sqrt_n, list_prime[i])); puts("");
	//exit(0);
	int lastE = 0;
	Int nextE = 1, las = sqrt_n, egol = sqrt_n * sqrt_n, fornext = nextE;
	for(int round = 0;; round++) {
		Int base = (Int) round * SIEVE_CHUNK + sqrt_n;
		fornext = nextE - base;

		for(register int i = 0; i < SIEVE_CHUNK; i++) {
			if (i >= fornext) {
				Int delta = i + base - las;
				Int Z = egol + (delta * delta) + (2 * delta * las);
				egol = Z;
				las = i + base;
				lastE = log2(Z - n);
				nextE = nextE * 9 / 5 + 1;
				fornext = nextE - base;
				temper[i] = lastE;
			}
			temper[i] = lastE;
		//	link[i].clear();
		}

		//cout<<temper[0]<<endl;
		rep(i, 0, list_size - 1) {
			int p = list_prime[i];
			float logp = log2((float)p);
			int res = list_res[i];
			int basemod = positive_mod(-base, p);
			int min_x = (basemod + res) % p;
			//printf("%d\n",min_x);
			while (min_x < SIEVE_CHUNK) {
				temper[min_x] -= logp;
				min_x += p;
			}
			if (p != 2) {
				res = p - res;
				min_x = (basemod + res) % p;
				while (min_x < SIEVE_CHUNK) {
					temper[min_x] -= logp;
					min_x += p;
				}
			}
			//cout<<"~"<<temper[0]<<" "<<p<<endl;
		}
		*nip = 0;
		rep(i, 0, SIEVE_CHUNK - 1) if (temper[i] < 1) {
			//printf("%d\n",i);
			nip[++*nip] = i;
		}
		rep(i, 0, SIEVE_CHUNK - 1) {
			temper[i] = -1;
		}
		rep(i, 1, *nip) {
			temper[nip[i]] = log2(y_f(nip[i], base, n));
			link[nip[i]].clear();
		}
		TT_clock -= clock();
		rep(i, 0, list_size - 1) {
			int p = list_prime[i];
			float logp = log2((float)p);
			int res = list_res[i];
			int min_x = positive_mod(res - base, p);
			//printf("%d\n",min_x);
			while (min_x < SIEVE_CHUNK) {
				if (temper[min_x] > 0) {
					temper[min_x] -= logp;
					link[min_x].push_back(i);
				}
				min_x += p;
			}
			if (p != 2) {
				res = p - res;
				min_x = positive_mod(res - base, p);
				while (min_x < SIEVE_CHUNK) {
					if (temper[min_x] > 0) {
						temper[min_x] -= logp;
						link[min_x].push_back(i);
					}
					min_x += p;
				}
			}
			//cout<<"~"<<temper[0]<<" "<<p<<endl;
		}
		TT_clock += clock();
		rep(i, 0, SIEVE_CHUNK - 1) if (temper[i] < 1 && temper[i] > 0) {
			//printf("%d\n", i);
			Ax.reset();
			dep(j, link[i].size() - 1, 0) {
			//	printf("%d ", link[i][j]);
				Ax[link[i][j]] = 1;
			}
			//puts("");
			id[A_cnt] = i + round * SIEVE_CHUNK;
			idAx.reset();
			idAx[A_cnt] = 1;
			T.reset();
			T[0] = 1;
			//cout<<(T & Ax)[0]<<endl;
			bool flag = true;
			for(int j = 0; j < PRIME_LIM; j++, T <<= 1) if ((T & Ax)[j]) {
				//if (A_cnt == 2) {
				//	printf("!A_cnt = %d, %d\n",A_cnt, i);
				//}
				if ((T & A[j])[j]) {
					Ax ^= A[j];
					idAx ^= idA[j];
				}
				else {
					A[j] = Ax;
					idA[j] = idAx;
					flag = false;
					break;
				}
			}
			if (flag) {
				//rep(j, 0, PRIME_LIM - 1) printf("%d ", (int)idAx[j]);
				//exit(0);

				//TT_clock += clock();
				if (check(sqrt_n, n, idAx)) {
				//TT_clock += clock();
					return;
				}
			}
			A_cnt++;
			printf("....%d/%d", A_cnt, PRIME_LIM);
			cout<<", number - sqrt{n} = "<<i + round * SIEVE_CHUNK<<endl;
		}

	}
}

int main() {
	prime_gen();
	Int inp;
	inp = "9999999999992989800000001205532601";
	//cout<<ln(inp)<<endl;
	//exit(0);
	int bee = clock();
	quadratic_sieve(inp);
	printf("time = %lf\n",(float)(clock()-bee)/1000000);
	printf("partial time = %'f\n", (float) TT_clock / 1000000);
	rep(i,0,ans.size()-1) cout<<ans[i]<<",";
	puts("");
	//rep(i, 0, 100) printf("%d ", prime[i]);
	//int z = quadratic_residue(7, 1000000007);
	//printf("%d", (long long) z * z % 1000000007);
}