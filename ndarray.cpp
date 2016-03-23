#include <stdio.h>
#include <stdlib.h>
#include "16807.h"

//Local Macros
#define POSITIVE 1
#define NEGATIVE -1

typedef char dtype;

struct ndarray{
    int dim;//dimension of the array
    int *shape;//shape of each dimension
    int *stride;//stride of
    dtype *head; //head of the array
};

void setmem(ndarray &arr, int dim, int *shape){
    int i=0;
    arr.dim=dim;
    arr.shape=(int *)malloc(sizeof(int)*arr.dim);
    arr.stride=(int *)malloc(sizeof(int)*(arr.dim+1));
    arr.stride[arr.dim]=1;
    for(i=0;i<arr.dim;i++){
        arr.shape[i]=shape[i];
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

void destroy(ndarray &arr){
    free(arr.shape);
    free(arr.stride);
    free(arr.head);
}

inline int offset(ndarray arr, int *index){
    int offcount=0,i;
    for(i=0;i<arr.dim;i++){
        offcount+=index[i]*arr.stride[i+1];
    }
    return offcount;
}

int getval(ndarray &arr, int *index){
    return arr.head[offset(arr,index)];
}

void setval(ndarray &arr, int *index, dtype val){
    arr.head[offset(arr,index)]=val;
}

inline int overflow(int index, int bound, int pn){
    if((pn == NEGATIVE)&&(index==0))
        return -1;
    else if((pn==POSITIVE) && (index==(bound-1)))
        return 1;
    return 0;
}

dtype rollval(ndarray arr, int rawind, int axis){
    int axisind, pn=POSITIVE;
    if(axis<0){
        pn=NEGATIVE;
        axis=-1-axis;
    }
    axisind=(rawind%arr.stride[axis])/arr.stride[axis+1];
    rawind+=pn*arr.stride[axis+1];
    rawind-=overflow(axisind, arr.shape[axis], pn)*arr.stride[axis];
    return arr.head[rawind];
}

int rollindex(ndarray arr, int rawind, int axis, int &dw){//copy of rollval
    int axisind, pn=POSITIVE;
    if(axis<0){
        pn=NEGATIVE;
        axis=-1-axis;
    }
    axisind=(rawind%arr.stride[axis])/arr.stride[axis+1];
    rawind+=pn*arr.stride[axis+1];
    dw=overflow(axisind, arr.shape[axis], pn);
    rawind-=dw*arr.stride[axis];
    return rawind;
}

void place(ndarray arr, double prob){
    for(int i=0;i<arr.stride[0];i++){
        if(myrand()<prob)
            arr.head[i]=1;
        else
            arr.head[i]=0;
    }
}

void printarr(ndarray &arr){
    int i,j;
    int ind[arr.dim+1]={0};
    if (arr.dim==1){
        printf("1D array:\n");
        for(i=0;i<arr.stride[0];i++)
            printf("%3d",arr.head[i]);
        putchar('\n');
        return;
    }
    else if(arr.dim==2){
        for(i=0;i<arr.shape[0];i++){
            for (j=0;j<arr.shape[1];j++){
                printf("%3d",arr.head[i*arr.stride[1]+j]);
            }
            putchar('\n');
        }
        return;
    }
    printf("Dimension:%d\tSize:%d\n",arr.dim,arr.stride[0]);
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
