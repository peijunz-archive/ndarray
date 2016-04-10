#ifndef NDARRAY_H
#define NDARRAY_H
#include <cstdio>
#include <iostream>
#include <cstdarg>
#include <initializer_list>
template <typename dtype>
/**
 * @brief The ndarray class
 * Defined the operations of a ndarray.
 */
class ndarray{
public:
    int dim;//dimension of the array
    int *shape;//shape of each dimension
    dtype *head; //head of the array
    ndarray():dim(0),shape(NULL),stride(NULL){}
    ndarray(int d, int width);
    ndarray(int d, int *sh);
    ndarray(std::initializer_list<int> l);
    void set(int d, int width);
    void set(int d, int *sh);
    ~ndarray();
    inline int size();
    void print();
    int rollindex(int rawind, int axis);
    dtype rollval(int rawind, int axis);
    inline dtype & operator[](int rawind);
    dtype & operator()(int *coo);
    dtype & operator()(int co0,...);
private:
    int *stride;//stride of every axis
    inline int overflow(int index, int bound, bool positive);
};

///Indexing elements
template <typename dtype>
/**
 * @brief ndarray<dtype>::operator []
 * @param offset    the raw index of the data
 * @return the element
 */
inline dtype & ndarray<dtype>::operator[](int offset){
    return head[offset];
}
template <typename dtype>
/**
 * @brief ndarray<dtype>::operator ()
 * @param coo,... the coordinates of an element
 * @return The indexed element
 *
 * Indexing by an index array
 */
dtype & ndarray<dtype>::operator()(int *coo){
    int offset=0;
    for(int i=0;i<dim;i++){
        offset+=stride[i+1]*coo[i];
    }
    return head[offset];
}
template <typename dtype>
/**
 * @brief ndarray<dtype>::operator ()
 * @param co0,... the index of an element
 * @return The indexed element
 *
 * Indexing by argument list
 */
dtype & ndarray<dtype>::operator()(int co0, ...){
    va_list l;
    va_start(l, co0);
    int offset=co0*stride[1];
    for(int i=1;i<dim;i++){
        offset+=stride[i+1]*va_arg(l,int);
    }
    va_end(l);
    return head[offset];
}

///Set the memory of ndarray
template <typename dtype>
/**
 * @brief ndarray<dtype>::set
 * @param d     dimension
 * @param w     width
 *
 * A square matrix
 */
void ndarray<dtype>::set(int d, int w){
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
template <typename dtype>
/**
 * @brief ndarray<dtype>::set
 * @param d     dimension
 * @param sh    shape array
 *
 * Set memory for an ndarray with a shape given by pointer/array
 */
void ndarray<dtype>::set(int d, int *sh){
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

///Initialize the ndarray by various method
template <typename dtype>
/**
 * @brief ndarray<dtype>::ndarray
 * @param l     initializer list for shape
 *
 * A ndarray with a shape given by an initializer list
 */
ndarray<dtype>::ndarray(std::initializer_list<int> l):dim(l.size()){
    shape=new int[dim];
    std::copy(std::begin(l), std::end(l), shape);
    stride=new int[dim+1];
    stride[dim]=1;
    for(int i=0;i<dim;i++){
        stride[dim-i-1]=stride[dim-i]*shape[dim-i-1];
    }
    head=new dtype[stride[0]];
}
template <typename dtype>
/**
 * @brief ndarray<dtype>::ndarray
 * @param d     dimension
 * @param w     width
 *
 * A square matrix
 */
ndarray<dtype>::ndarray(int d, int w){
    set(d, w);
}
template <typename dtype>
/**
 * @brief ndarray<dtype>::ndarray
 * @param d     dimension
 * @param sh    shape array
 *
 * A ndarray with a shape given by pointer/array
 */
ndarray<dtype>::ndarray(int d, int *sh){
    set(d, sh);
}

template <typename dtype>
/**
 * @brief ndarray<dtype>::~ndarray
 * Destructor
 */
ndarray<dtype>::~ndarray(){
    delete [] shape;
    delete [] stride;
    delete [] head;
}

template <typename dtype>
/**
 * @brief ndarray<dtype>::size
 * @return The size of data
 */
inline int ndarray<dtype>::size(){
    return stride[0];
}

template <typename dtype>
inline int ndarray<dtype>::overflow(int index, int bound, bool positive){
    if((!positive)&&(index==0))
        return -1;
    else if(positive && (index==(bound-1)))
        return 1;
    return 0;
}
template <typename dtype>
int ndarray<dtype>::rollindex(int rawind, int axis){//copy of rollval
    int axisind;
    bool pn;
    if(axis<0){
        pn=false;
        axis=-1-axis;
    }
    else pn=true;
    axisind=(rawind%stride[axis])/stride[axis+1];
    rawind+=pn*stride[axis+1];
    rawind-=overflow(axisind, shape[axis], pn)*stride[axis];
    return rawind;//Compare new raw index will indicate overflow or not
}
template <typename dtype>
dtype ndarray<dtype>::rollval(int rawind, int axis){
    return head[rollindex(rawind,axis)];
}

template <typename dtype>
/**
 * @brief ndarray<dtype>::print print the matrix
 *
 * For 1D array, print horizontally
 * For 2D array, printed as an matrix
 * For nD array(n>2), printed all indices and corresponding values.
 */
void ndarray<dtype>::print(){
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

template <typename T>
/**
 * @brief homogenize    set all the values in the array to a init value
 * @param a         array
 * @param initval   initial value
 */
void homogenize(ndarray<T> &a, T initval=0){
    for(int i=0;i<a.size();i++)
        a[i]=initval;
    return;
}

#endif //NDARRAY_H
