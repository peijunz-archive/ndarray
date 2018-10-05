#include <iostream>
#include "ndarray.h"
#include "nditer.h"
using namespace std;

int main()
{
    int m=10000, n=1024;
    matrix<int> arr{m, n};
    arr = 1;
    for (int i=1; i<m; i++){
        for (int j=1; j<n; j++){
//             arr[i*arr.stride(0)+j] += arr[(i-1)*arr.stride(0)+(j-1)] + arr[(i-1)*arr.stride(0)+j] + arr[i*arr.stride(0)+(j-1)];
//             arr[i*n+j] += arr[(i-1)*n+(j-1)] + arr[(i-1)*n+j] + arr[i*n+(j-1)];
//             arr[i*1024+j] += arr[(i-1)*1024+(j-1)] + arr[(i-1)*1024+j] + arr[i*1024+(j-1)];
            arr(i, j) += arr((i-1), (j-1)) + arr((i-1), j) + arr(i, j-1);
        }
    }
    cout << arr(m-1, n-1)<<endl;
    cout<<arr;
    return 0;
}
