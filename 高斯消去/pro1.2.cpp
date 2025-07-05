#include<iostream>
#include<Windows.h>
using namespace std;

int main() {
    int n;
    int k;
    cout << "��������Թ�ģ��" << endl;
    cin >> k;
    cout << "����������ģ��" << endl;
    cin >> n;

    // ����n*n����
    int** b = new int* [n];
    for (int i = 0; i < n; i++) {
        b[i] = new int[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i][j] = i + j;
        }
    }

    // ��������
    int* a = new int[n];
    // �����Ž��������
    int* sum = new int[n];
    for (int i = 0; i < n; i++) {
        a[i] = i;
    }

    long long head, tail, freq; // timers
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);


    // cache�Ż��㷨�����������ڻ�
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
    cout << "cache�Ż��㷨��ʱ��Ϊ��" << (tail - head) * 1000.0 / freq << "ms" << endl;
    cout << "cache�Ż��㷨ƽ��ʱ��Ϊ:" << (tail - head) * 1000.0 / freq / k << "ms" << endl;

    // �ͷ��ڴ�
    delete[]a;
    for (int i = 0; i < n; i++) {
        delete[]b[i];
    }
    delete[]b;
    delete[]sum;

    return 0;
}