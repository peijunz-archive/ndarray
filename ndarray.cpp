#include <stdio.h>
#include <stdlib.h>
#include "16807.h"
#define NORMAL 0
#define SHIFT 1
#define SINGLE 2
typedef char dtype;
struct ndarray{
    int dim;//dimension of the array
    int flag;
    int *shape;//shape of each dimension
    int *stride;//stride of
    int *shift;//shift of a specific axis
    dtype* head; //head of the array
};
//inline int size(ndarray arr){
//    return arr.stride[0];
//}
void printarr(ndarray arr){
    int i,j;
    int ind[arr.dim+1]={0};
    printf("Dimension:%d\tSize:%d\tFlag:%d\nElements:\n",arr.dim,arr.stride[0],arr.flag);
    printf("Stride\t");
    for(i=0;i<arr.dim;i++)
        printf("%d\t",arr.stride[i+1]);
    printf("\nShape\t");
    for(i=0;i<arr.dim;i++)
        printf("%d\t",arr.shape[i]);
    printf("\nRaw\t");
    for(i=0;i<arr.dim;i++)
        printf("Axis %d\t",i);
    printf("Value\n");
    for(i=0;i<arr.stride[0];i++){
        printf("%d\t",i);
        for(j=arr.dim-1;ind[j+1]==arr.shape[j];j--){
            ind[j+1]=0;
            ind[j]+=1;
        }
        for(j=0;j<arr.dim;j++)
            printf("%d\t",ind[j+1]);
        printf("%d\n",arr.head[i]);
        ind[arr.dim]+=1;
    }
}

void setmem(ndarray &arr, int dim, int *shape){
    int i=0;
    arr.dim=dim;
    arr.shape=(int *)malloc(sizeof(int)*arr.dim);
    arr.stride=(int *)malloc(sizeof(int)*(arr.dim+1));
    arr.shift=(int *)malloc(sizeof(int)*arr.dim);
    arr.stride[arr.dim]=1;
    for(i=0;i<arr.dim;i++){
        arr.shape[i]=shape[i];
        arr.shift[i]=0;
        arr.stride[dim-i-1]=arr.stride[dim-i]*shape[dim-i-1];
    }
    arr.head=(dtype *)malloc(sizeof(dtype)*arr.stride[0]);
    return;
}
void init(ndarray &arr, dtype initval=0){
    for(int i=0;i<arr.stride[0];i++)
        arr.head[i]=initval;
    return;
}
void init(ndarray &arr, dtype *initarr){
    for(int i=0;i<arr.stride[0];i++)
        arr.head[i]=initarr[i];
    return;
}

void destroy(ndarray &arr){
    free(arr.shape);
    free(arr.stride);
    free(arr.shift);
    free(arr.head);
}

inline int modnav(int a, int b){
    //Naive modulus function
    if (a>b-1)
        return a-b;
    else if(a<0)
        return a+b;
    return a;
}
inline int offset(ndarray arr, int *index, int flag=0){
    int offcount=0,i;
    //for general condition, offset=stride*((index+shift)mod(shape))
    if (flag==0){
        for(i=0;i<arr.dim;i++){
            offcount+=index[i]*arr.stride[i+1];
        }
    }
    else if(flag==SHIFT){
        for(i=0;i<arr.dim;i++){
            offcount+=modnav((index[i]+arr.shift[i]),arr.shape[i])*arr.stride[i+1];
        }
    }
    return offcount;
}

int getval(ndarray &arr, int *index, int flag=0){
    return arr.head[offset(arr,index,flag)];
}
int setval(ndarray &arr, int *index, dtype val, int flag=0){
    arr.head[offset(arr,index,flag)]=val;
}
inline int overflow(int index, int bound, bool pn){
    if((!pn)&&(index==0))
        return -1;
    else if(pn && (index==(bound-1)))
        return 1;
    return 0;
}
dtype rollval(ndarray arr, int rawind, int axis, bool pn=true){
    int axisind=(rawind%arr.stride[axis])/arr.stride[axis+1];
    rawind+=arr.stride[axis+1];
    rawind-=overflow(axisind, arr.shape[axis], pn)*arr.stride[axis];
    return arr.head[rawind];
}

void place(ndarray arr, double prob=0.5){
    for(int i=0;i<arr.stride[0];i++){
        if(myrand()<prob)
            arr.head[i]=1;
        else
            arr.head[i]=0;
    }
}
