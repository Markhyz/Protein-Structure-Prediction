#include <bits/stdc++.h>

using namespace std;

int main() {
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> urd(0, 1);
    cout << fixed;
    for (int i = 0; i < 1e9; ++i) {
        double x = i / 1e15;
        double y = x * 1e15;
        if (fabs(i - y) > 1e-7) {
            cout.precision(8);
            cout << i << " " << x << " " << x * 1e15 << " " << y << endl;
            exit(-1);
        }
    }

    cout << "ok" << endl;

    return 0;
}
