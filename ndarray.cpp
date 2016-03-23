#include <stdio.h>
#include <stdlib.h>
#include "16807.h"

//Local Macros
#define POSITIVE 1
#define NEGATIVE -1

typedef char dtype;

class ndarray{
public:
    int dim;//dimension of the array
    int *shape;//shape of each dimension
    int *stride;//stride of
    dtype *head; //head of the array
    ndarray():dim(0),shape(NULL),stride(NULL){}
    ndarray(int d, int width);
    ndarray(int d, int *sh);
    void reset(int d, int width);
    void reset(int d, int *sh);
    ~ndarray();
    inline int size();
    void print();
    int rollindex(int rawind, int axis);
    dtype rollval(int rawind, int axis);
private:
    inline int overflow(int index, int bound, int pn);
};

void ndarray::reset(int d, int w){
    if (dim>0){
        delete [] shape;
        delete [] stride;
        delete [] head;
    }
    dim=d;
    shape=new int[d];
    stride=new int[d+1];
    stride[d]=1;
    for(int i=0;i<d;i++){
        shape[i]=w;
        stride[d-i-1]=stride[d-i]*w;
    }
    head=new dtype[stride[0]];
}
void ndarray::reset(int d, int *sh){
    if (dim>0){
        delete [] shape;
        delete [] stride;
        delete [] head;
    }
    dim=d;
    shape=new int[d];
    stride=new int[d+1];
    stride[d]=1;
    for(int i=0;i<d;i++){
        shape[i]=sh[i];
        stride[d-i-1]=stride[d-i]*sh[d-i-1];
    }
    head=new dtype[stride[0]];
}

ndarray::ndarray(int d, int *sh):dim(d){
    shape=new int[d];
    stride=new int[d+1];
    stride[d]=1;
    for(int i=0;i<d;i++){
        shape[i]=sh[i];
        stride[d-i-1]=stride[d-i]*sh[d-i-1];
    }
    head=new dtype[stride[0]];
}
ndarray::ndarray(int d, int w):dim(d){
    shape=new int[d];
    stride=new int[d+1];
    stride[d]=1;
    for(int i=0;i<d;i++){
        shape[i]=w;
        stride[d-i-1]=stride[d-i]*w;
    }
    head=new dtype[stride[0]];
}
ndarray::~ndarray(){
    delete [] shape;
    delete [] stride;
    delete [] head;
}

inline int ndarray::size(){
    return stride[0];
}

inline int ndarray::overflow(int index, int bound, int pn){
    if((pn == NEGATIVE)&&(index==0))
        return -1;
    else if((pn==POSITIVE) && (index==(bound-1)))
        return 1;
    return 0;
}

int ndarray::rollindex(int rawind, int axis){//copy of rollval
    int axisind, pn;
    if(axis<0){
        pn=NEGATIVE;
        axis=-1-axis;
    }
    else pn=POSITIVE;
    axisind=(rawind%stride[axis])/stride[axis+1];
    rawind+=pn*stride[axis+1];
    rawind-=overflow(axisind, shape[axis], pn)*stride[axis];
    return rawind;//Compare new raw index will indicate overflow or not
}

dtype ndarray::rollval(int rawind, int axis){
    return head[rollindex(rawind,axis)];
}

void ndarray::print(){
    int i,j;
    int ind[dim+1]={0};
    if(dim==0){
        printf("NO ARRAY!\n");
    }
    else if (dim==1){
        printf("1D array:\n");
        for(i=0;i<stride[0];i++)
            printf("%3d",head[i]);
        putchar('\n');
        return;
    }
    else if(dim==2){
        for(i=0;i<shape[0];i++){
            for (j=0;j<shape[1];j++){
                printf("%3d",head[i*stride[1]+j]);
            }
            putchar('\n');
        }
        return;
    }
    printf("Dimension:%d\tSize:%d\n",dim,stride[0]);
    printf("Stride\t");
    for(i=0;i<dim;i++)
        printf("%d\t",stride[i+1]);
    printf("\nShape\t");
    for(i=0;i<dim;i++)
        printf("%d\t",shape[i]);
    printf("\nRaw\t");
    for(i=0;i<dim;i++)
        printf("Axis %d\t",i);
    printf("Value\n");
    for(i=0;i<stride[0];i++){
        printf("%d\t",i);
        for(j=dim-1;ind[j+1]==shape[j];j--){
            ind[j+1]=0;
            ind[j]+=1;
        }
        for(j=0;j<dim;j++)
            printf("%d\t",ind[j+1]);
        printf("%d\n",head[i]);
        ind[dim]+=1;
    }
}

void randomize(ndarray &a, double prob){
    for(int i=0;i<a.stride[0];i++){
        if(myrand()<prob)
            a.head[i]=1;
        else
            a.head[i]=0;
    }
}

void homogenize(ndarray &a, dtype initval=0){
    for(int i=0;i<a.stride[0];i++)
        a.head[i]=initval;
    return;
}

//inline int offset(ndarray arr, int *index){
//    int offcount=0,i;
//    for(i=0;i<arr.dim;i++){
//        offcount+=index[i]*arr.stride[i+1];
//    }
//    return offcount;
//}

//int getval(ndarray &arr, int *index){
//    return arr.head[offset(arr,index)];
//}

//void setval(ndarray &arr, int *index, dtype val){
//    arr.head[offset(arr,index)]=val;
//}
//void setmem(ndarray &arr, int dim, int *shape){
//    int i=0;
//    arr.dim=dim;
//    arr.shape=(int *)malloc(sizeof(int)*arr.dim);
//    arr.stride=(int *)malloc(sizeof(int)*(arr.dim+1));
//    arr.stride[arr.dim]=1;
//    for(i=0;i<arr.dim;i++){
//        arr.shape[i]=shape[i];
//        arr.stride[dim-i-1]=arr.stride[dim-i]*shape[dim-i-1];
//    }
//    arr.head=(dtype *)malloc(sizeof(dtype)*arr.stride[0]);
//    return;
//}

//void destroy(ndarray &arr){
//    free(arr.shape);
//    free(arr.stride);
//    free(arr.head);
//}

int prod(int dim, int *shape){
    int size=1;
    for(int i=0;i<dim;i++){
        size*=shape[i];
    }
    return size;
}
int intpow(int x, int y){
    int res=1;
    for(int i=0;i<y;i++){
        res*=x;
    }
    return res;
}
