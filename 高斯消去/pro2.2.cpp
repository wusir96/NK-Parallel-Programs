#include <iostream>
#include <Windows.h>
using namespace std;



// 适合超标量架构的指令级并行算法
int parallelSum(int* a, int n) {
    int sum1 = 0, sum2 = 0;
    for (int i = 0; i < n; i += 2) {
        sum1 += a[i];
        if (i + 1 < n) {
            sum2 += a[i + 1];
        }
    }
    return sum1 + sum2;
}

int main() {
    cout << "请输入测试规模：" << endl;
    int k;
    cin >> k;
    cout << "请输入数组规模：" << endl;
    int n;
    cin >> n;

    // 定义累加的数组
    int* a = new int[n];
    for (int i = 0; i < n; i++) {
        a[i] = i;
    }

    // 定义累加结果
    int sum;
    long long head, tail, freq; // timers

    // 获取计时器频率
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);


    // 适合超标量架构的指令级并行算法计时
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    for (int cnt = 0; cnt < k; cnt++) {
        sum = parallelSum(a, n);
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "指令级并行算法总时间为：" << (tail - head) * 1000.0 / freq << "ms" << endl;
    cout << "指令级并行算法平均时间为：" << (tail - head) * 1000.0 / freq / k << "ms" << endl;

    // 释放内存
    delete[] a;
    return 0;
}