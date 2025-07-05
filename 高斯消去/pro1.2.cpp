#include<iostream>
#include<Windows.h>
using namespace std;

int main() {
    int n;
    int k;
    cout << "请输入测试规模：" << endl;
    cin >> k;
    cout << "请输入矩阵规模：" << endl;
    cin >> n;

    // 定义n*n矩阵
    int** b = new int* [n];
    for (int i = 0; i < n; i++) {
        b[i] = new int[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i][j] = i + j;
        }
    }

    // 定义向量
    int* a = new int[n];
    // 定义存放结果的向量
    int* sum = new int[n];
    for (int i = 0; i < n; i++) {
        a[i] = i;
    }

    long long head, tail, freq; // timers
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);


    // cache优化算法矩阵与向量内积
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    for (int i = 0; i < k; i++) {
        for (int cnt = 0; cnt < n; cnt++) {
            sum[cnt] = 0;
        }
        for (int col = 0; col < n; col++) {
            for (int row = 0; row < n; row++) {
                sum[row] += a[col] * b[row][col];
            }
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "cache优化算法总时间为：" << (tail - head) * 1000.0 / freq << "ms" << endl;
    cout << "cache优化算法平均时间为:" << (tail - head) * 1000.0 / freq / k << "ms" << endl;

    // 释放内存
    delete[]a;
    for (int i = 0; i < n; i++) {
        delete[]b[i];
    }
    delete[]b;
    delete[]sum;

    return 0;
}