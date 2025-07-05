#include <iostream>
#include <Windows.h>
using namespace std;

// ����ۼӵ�ƽ���㷨
int trivialSum(int* a, int n) {
    int sum = 0;
    for (int i = 0; i < n; i++) {
        sum += a[i];
    }
    return sum;
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

    // ����ۼӵ�ƽ���㷨��ʱ
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    for (int cnt = 0; cnt < k; cnt++) {
        sum = trivialSum(a, n);
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "����ۼӵ�ƽ���㷨��ʱ��Ϊ��" << (tail - head) * 1000.0 / freq << "ms" << endl;
    cout << "����ۼӵ�ƽ���㷨ƽ��ʱ��Ϊ��" << (tail - head) * 1000.0 / freq / k << "ms" << endl;

    
    // �ͷ��ڴ�
    delete[] a;
    return 0;
}