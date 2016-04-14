#ifndef NDARRAY_H
#define NDARRAY_H
#include <iostream>
#include <cstdarg>
#include <initializer_list>
using namespace std;
/**
 * @file ndarray.h
 * @author zpj
 * @brief The n-Dimensional array template
 * @bug No
 */
template <typename dtype>
/**
 * @brief The n-Dimensional array template class
 */
class ndarray{
public:
    int dim;            ///< dimension of the array
    int *shape;         ///< shape of each dimension
    dtype *head;        ///< head of the array data
    /// Initialization by setting all data to zero without allocating memory
    ndarray():dim(0),shape(NULL),head(NULL),stride(NULL){}
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
    void homogenize(dtype initval);
private:
    int *stride;        ///< stride of every axis
    inline int overflow(int index, int bound, bool positive);
};

template <typename dtype>
/**
 * @brief Indexing by raw index
 * @param offset    the raw index of the data
 * @return the element
 */
inline dtype & ndarray<dtype>::operator[](int offset){
    return head[offset];
}
template <typename dtype>
/**
 * @brief Indexing by an index array
 * @param coo,... the coordinates of an element
 * @return The indexed element
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
 * @brief Indexing by the argument list of operator()
 * @param co0,... the index of an element
 * @return The indexed element
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

template <typename dtype>
/**
 * @brief Set memory for an a square matrix
 * @param d     dimension
 * @param w     width
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
 * @brief Set memory for an ndarray with a shape given by pointer/array
 * @param d     dimension
 * @param sh    shape array
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

template <typename dtype>
/**
 * @brief Initialize an ndarray with a shape given by an initializer list
 * @param l     initializer list for shape
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
 * @brief Initialize an ndarray for a square matrix
 * @param d     dimension
 * @param w     width
 */
ndarray<dtype>::ndarray(int d, int w){
    set(d, w);
}
template <typename dtype>
/**
 * @brief A ndarray with a shape given by pointer/array
 * @param d     dimension
 * @param sh    shape array
 */
ndarray<dtype>::ndarray(int d, int *sh){
    set(d, sh);
}

template <typename dtype>
/**
 * @brief Free all allocated memory.
 */
ndarray<dtype>::~ndarray(){
    delete [] shape;
    delete [] stride;
    delete [] head;
}

template <typename dtype>
/**
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
/**
 * @brief Roll the index in a given axis with periodical boundary condition
 * @param rawind raw index
 * @param axis
 * @return Raw index after rolling
 */
int ndarray<dtype>::rollindex(int rawind, int axis){
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
    return rawind;
}
template <typename dtype>
/**
 * @brief Roll the index in a given axis with periodical boundary condition
 * @param rawind raw index
 * @param axis
 * @return Corresponding value after roll
 */
dtype ndarray<dtype>::rollval(int rawind, int axis){
    return head[rollindex(rawind,axis)];
}

template <typename dtype>
/**
 * @brief print the matrix/ndarray
 *
 * + 1D array, printed horizontally
 * + 2D array, printed as an matrix
 * + 3+D array, printed all indices and corresponding values.
 */
void ndarray<dtype>::print(){
    int i,j;
    int ind[dim+1]={0};
    if(dim==0){
        cout<<"NO ARRAY!"<<endl;
    }
    else if (dim==1){
        cout<<"1D array:"<<endl;
        for(i=0;i<stride[0];i++)
            cout<<head[i]<<'\t';
        cout<<endl;
        return;
    }
    else if(dim==2){
        for(i=0;i<shape[0];i++){
            for (j=0;j<shape[1];j++){
                cout<<head[i*stride[1]+j]<<'\t';
            }
            cout<<endl;
        }
        return;
    }
    cout<<"Dimension:"<<dim<<'\t'<<"Size:"<<stride[0]<<'\n';
    cout<<"Stride\t";
    for(i=0;i<dim;i++)
        cout<<stride[i+1]<<'\t';
    cout<<"\nShape\t";
    for(i=0;i<dim;i++)
        cout<<shape[i]<<"\t";
    cout<<"\nRaw\t";
    for(i=0;i<dim;i++)
        cout<<"Axis "<<i<<"\t";
    cout<<"Value\n";
    for(i=0;i<stride[0];i++){
        cout<<i<<"\t";
        for(j=dim-1;ind[j+1]==shape[j];j--){
            ind[j+1]=0;
            ind[j]+=1;
        }
        for(j=0;j<dim;j++)
            cout<<ind[j+1]<<"\t";
        cout<<head[i]<<endl;
        ind[dim]+=1;
    }
}

template <typename dtype>
/**
 * @brief set all the values in the array to a init value
 * @param initval   initial value
 */
void ndarray<dtype>::homogenize(dtype initval){
    for(int i=0;i<size();i++)
        head[i]=initval;
    return;
}


#endif //NDARRAY_H
