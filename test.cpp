#include <iostream>
#include "ndarray.h"
#include "nditer.h"
using namespace std;

int main()
{
    ndarray<int, 2> arr{20,10};
    auto coo=nditer<2>{20,10};
    for(int i=0;i<arr.size();i++){
        arr[i]=i;
    }
    for(const auto k: arr){
        cout<<k<<endl;
    }
    do{
        coo.print();
        cout<<arr(coo())<<endl;
    }while(++coo);
    cout<<sizeof (arr)<<endl;
    return 0;
}
