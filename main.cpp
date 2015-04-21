#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstring>
#include <string>
#include <algorithm>
#include <pthread.h>
#define CPS CLOCKS_PER_SEC

class bits {
public:
  unsigned long a[4];
  void clear() {
    a[0] = a[1] = a[2] = a[3] = 0;
  }
  void set(int x) {
    a[x / 64] |= 1ul << (x % 64);
  }
  void reset(int x) {
    a[x / 64] &= (~0ul) ^ (1 << (x % 64));
  }
  friend bits& operator <<= (bits &a, const int b) {
    for (int i = 3; i >= 0; i--) {
      a.a[i] = (a.a[i] << b) | (i ? (a.a[i - 1] >> (64 - b)) : 0);
    }
    return a;
  }
  friend bool operator < (const bits &a, const bits &b) {
    int i;
    for (i = 0; i < 3 && a.a[i] == b.a[i]; i++);
    return a.a[i] < b.a[i];
  }
};

char dna[1048576][128], buf[128];
int lst[1048576 * 128];
int n, m, s, threads = 8;
bits hash[1048576][128];

void load_file(std::string filename) {
  FILE *file = fopen(filename.c_str(), "r");

  while (fgets(buf, 128, file) != NULL)
    assert(strlen(fgets(dna[n++], 128, file)) == 101);

  fclose(file);
}

void read_all_dnas() {
  clock_t st = clock();

  load_file("solexa_100_170_1.fa");
  load_file("solexa_100_170_2.fa");
  printf("%d dns loaded\n", n);

  clock_t ed = clock();
  printf("--->load file time: %.6lfs\n", (double)(ed - st) / CPS);
}

bool cmp_lst(int x, int y) {
  bits *a = &hash[x / 128][x & 127], *b = &hash[y / 128][y & 127];
  return *a < *b;
}

void *preprocess_thread(void *id) {
  for (long i = (long)id; i < n; i += threads) {
    for (int j = 0; j < m; j++) {
      hash[i][0] <<= 2;
      if (dna[i][j] == 'A' || dna[i][j] == 'C')
        hash[i][0].set(1);
      if (dna[i][j] == 'A' || dna[i][j] == 'G')
        hash[i][0].set(0);
    }
    for (int j = 1; j < 101 - m; j++) {
      hash[i][j] = hash[i][j - 1];
      hash[i][j] <<= 2;
      hash[i][j].reset(2 * m);
      hash[i][j].reset(2 * m + 1);
      if (dna[i][j + m - 1] == 'A' || dna[i][j + m - 1] == 'C')
        hash[i][j].set(1);
      if (dna[i][j + m - 1] == 'A' || dna[i][j + m - 1] == 'G')
        hash[i][j].set(0);
    }
  }
  return NULL;
}

void preprocess() {
  printf("Input K: ");
  scanf("%d", &m);

  clock_t st = clock();

  pthread_t *pt = (pthread_t *)malloc(threads * sizeof(pthread_t));
  for (long i = 0; i < threads; i++) pthread_create(&pt[i], NULL, preprocess_thread, (void *)i);
  for (long i = 0; i < threads; i++) pthread_join(pt[i], NULL);

  clock_t ed = clock();
  float tt = (double)(ed - st) / CPS / threads;

  st = clock();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 101 - m; j++)
      lst[s++] = i * 128 | j;

  std::sort(lst, lst + s, cmp_lst);
  ed = clock();
  printf("%d k-murs\n", s);
  printf("--->pre-process time: %.6lfs\n", (double)(ed - st) / CPS + tt);
}

bool cmp_l(const int x, const bits &y) {
  bits *a = &hash[x / 128][x & 127];
  return *a < y;
}

bool cmp_u(const bits &x, const int y) {
  bits *b = &hash[y / 128][y & 127];
  return x < *b;

}

void solve() {
  printf("Input a pattern: ");
  scanf("%s", buf);
  assert(strlen(buf) == m);

  float st = clock();

  bits opt;
  opt.clear();
  for (int i = 0; i < m; i++) {
    opt <<= 2;
    if (buf[i] == 'A' || buf[i] == 'C')
      opt.set(1);
    if (buf[i] == 'A' || buf[i] == 'G')
      opt.set(0);
  }

  long x = std::lower_bound(lst, lst + s, opt, cmp_l) - lst;
  long y = std::upper_bound(lst, lst + s, opt, cmp_u) - lst;
  printf("%ld k-mur found.\n", y - x);

  float ed = clock();
  printf("--->solve 1 query time: %.6lfs\n", (double)(ed - st) / CPS);
  for (int i = x; i < y; i++) {
    printf("string %d, pos %d. \n", lst[i] / 128 + 1, (lst[i] & 127) + 1);
  }
}

int main(int argc, char **argv) {
  if (argc == 2)
    threads = atoi(argv[1]);

  read_all_dnas();

  preprocess();

  while (1)
    solve();
  return 0;
}
