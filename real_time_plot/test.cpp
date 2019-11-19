#include <unistd.h>
#include <cstdio>
#include <random>

void g(double x, double y) {
  printf("nani %.5f\n", x * y);
}

void g(int x) {
  printf("yoo %d\n", x);
}

template<typename ...ArgsT>
void f(ArgsT... args) {
  g(args...);
}

int main() {

  FILE* fp = popen("python test.py", "w");
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0,50);

  for (int i = 0; i < 1000000; ++i) {
    int x, y; 
    x = i; y = distribution(generator);
    fprintf(fp, "%d %d\n", x, y);
    fflush(fp);
    usleep(100000);
  }

  return 0;
}
