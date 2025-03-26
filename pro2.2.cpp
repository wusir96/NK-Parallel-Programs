#include <iostream>
#include <Windows.h>
using namespace std;



// �ʺϳ������ܹ���ָ������㷨
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
    cout << "��������Թ�ģ��" << endl;
    int k;
    cin >> k;
    cout << "�����������ģ��" << endl;
    int n;
    cin >> n;

    // �����ۼӵ�����
    int* a = new int[n];
    for (int i = 0; i < n; i++) {
        a[i] = i;
    }

    // �����ۼӽ��
    int sum;
    long long head, tail, freq; // timers

    // ��ȡ��ʱ��Ƶ��
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);


    // �ʺϳ������ܹ���ָ������㷨��ʱ
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    for (int cnt = 0; cnt < k; cnt++) {
        sum = parallelSum(a, n);
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "ָ������㷨��ʱ��Ϊ��" << (tail - head) * 1000.0 / freq << "ms" << endl;
    cout << "ָ������㷨ƽ��ʱ��Ϊ��" << (tail - head) * 1000.0 / freq / k << "ms" << endl;

    // �ͷ��ڴ�
    delete[] a;
    return 0;
}