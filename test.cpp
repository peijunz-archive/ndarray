#include <iostream>
#include "ndarray.h"
#include "nditer.h"
using namespace std;

int main()
{
    ndarray<int, 2> arr{20,10};
    auto coo=arr.index();
    for(int i=0;i<arr.size();i++){
        arr[i]=i;
    }
    do{
        coo.print();
        cout<<arr(coo())<<endl;
    }while(++coo);
    return 0;
}
