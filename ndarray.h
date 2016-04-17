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
/// The n-Dimensional array template class
class ndarray{
private:
    int dim;            ///< dimension of the array
    int *shape;         ///< shape of each dimension
    int *stride;        ///< stride of every axis
public:
    dtype *head;        ///< head of the array data

    //Constructors and Assignments
    ndarray():dim(0),shape(nullptr),stride(&dim),head(nullptr){}
    ndarray(int d, int width);
    ndarray(int d, int *sh);
    ndarray(initializer_list<int> l);
    ndarray(const ndarray<dtype>&);
    ndarray(ndarray<dtype> && x);
    ~ndarray<dtype>();
    ndarray<dtype> & operator=(const ndarray<dtype> &);
    ndarray<dtype> & operator=(ndarray<dtype> && x);
    ndarray<dtype> & operator=(dtype val);

    //Read private data
    int size() {return stride[0];}      ///< The size of data as an mask of stride[0]
    int _shape(int i){return shape[i];}    ///< Get shape
    int _dim(){return dim;}                ///< Get dimension
    int _stride(int i){return stride[i+1];}///< Get stride

    //Indexing Data
    inline dtype & operator [](int rawind);
    dtype & operator ()(int *coo);
    dtype & operator ()(int co0,...);
    int rollindex(int rawind, int axis);
    //dtype rollval(int rawind, int axis);
    //Other
    void print();
private:
    inline int overflow(int index, int bound, bool positive);
};
//Constructor
//template <typename dtype>
///**
// * @brief Default constructor
// */
//ndarray<dtype>::ndarray():
//    dim(0),shape(nullptr),head(nullptr),stride(new int[1])
//{
//    stride=&dim;
//}
template <typename dtype>
/**
 * @brief Initialize an ndarray for a square matrix
 * @param d     dimension
 * @param w     width
 */
ndarray<dtype>::ndarray(int d, int w):
    dim(d),shape(new int[d]),stride(new int[d+1]){
    stride[d]=1;
    for(int i=0;i<d;i++){
        shape[i]=w;
        stride[d-i-1]=stride[d-i]*w;
    }
    head=new dtype[stride[0]];
}
template <typename dtype>
/**
 * @brief A ndarray with a shape given by pointer/array
 * @param d     dimension
 * @param sh    shape array
 *
 * For dynamic construction
 */
ndarray<dtype>::ndarray(int d, int *sh):
    dim(d),shape(new int[d]),stride(new int[d+1]){
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
 *
 * For static construction
 */
ndarray<dtype>::ndarray(initializer_list<int> l):dim(l.size()){
    shape=new int[dim];
    copy(begin(l), end(l), shape);
    stride=new int[dim+1];
    stride[dim]=1;
    for(int i=0;i<dim;i++){
        stride[dim-i-1]=stride[dim-i]*shape[dim-i-1];
    }
    head=new dtype[stride[0]];
}
//Destructor
template <typename dtype>
/// Free all allocated memory.
ndarray<dtype>::~ndarray(){
    if(dim){
        delete [] shape;
        delete [] stride;
        delete [] head;
    }
}
//Deep Copy. not copy data should be aliased as "like"
template <typename dtype>
/**
 * @brief Copy constructor copying the shape only.
 * @param x     The source of ndarray
 */
ndarray<dtype>::ndarray(const ndarray<dtype> & x):
    dim(x.dim),shape(new int[x.dim]),\
    stride(new int[x.dim+1]),\
    head(new dtype[x.stride[0]])
{
    stride[dim]=1;
    for(int i=0;i<dim;i++){
        shape[i]=x.shape[i];
        stride[i]=x.stride[i];
    }
}

template <typename dtype>
/**
 * @brief Copy Assignment copying all data
 * @param x     The source of ndarray
 */
ndarray<dtype> & ndarray<dtype>:: operator= (const ndarray<dtype>& x){
    if(x.stride[0]==0){
        exit(0);
        cout<<"Error: Copy an empty ndarray! Use move constructor instead.";
    }
    if(stride[0]!=x.stride[0]){
        delete [] head;
        head=new dtype[x.stride[0]];
    }
    if(dim==x.dim){
        bool unequal=false;
        for(int i=0;i<dim;i++){
            if(shape[i]!=x.shape[i]){
                unequal=true;
                shape[i]=x.shape[i];
            }
        }
        if(unequal){
            for(int i=0;i<dim;i++){
                stride[i]=x.stride[i];
            }
        }
    }
    else{
        if(dim){//In case that in dim=0 case, delete [] &dim
            delete [] shape;
            delete [] stride;
        }
        dim=x.dim;
        shape=new int[x.dim];
        stride=new int[x.dim+1];
        stride[dim]=1;
        for(int i=0;i<dim;i++){
            stride[i]=x.stride[i];
            shape[i]=x.shape[i];
        }
    }
    for(int i=0;i<stride[0];i++){
        head[i]=x.head[i];
    }
    return *this;
}
//Move
template <typename dtype>
/// Move constructor
ndarray<dtype>::ndarray(ndarray<dtype> && x):
    dim(x.dim),shape(x.shape),head(x.head),stride(x.stride){
    x.shape=nullptr;
    x.head=nullptr;
    x.stride=nullptr;
}
template <typename dtype>
/// Move assignment
ndarray<dtype> & ndarray<dtype>:: operator= (ndarray<dtype>&& x){
    dim=x.dim;
    shape=x.shape;
    head=x.head;
    stride=x.stride;
    x.shape=nullptr;
    x.head=nullptr;
    x.stride=nullptr;
    return *this;
}
//Indexing
template <typename dtype>
/**
 * @brief Indexing by *raw* index
 * @param offset    the raw index of the data
 * @return the element
 */
inline dtype & ndarray<dtype>:: operator [](int offset){
    return head[offset];
}
template <typename dtype>
/**
 * @brief Indexing by an index array
 * @param coo,... the coordinates of an element
 *
 * Especially useful for array with known dimension.
 * @return The indexed element
 */
dtype & ndarray<dtype>:: operator ()(int *coo){
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
 *
 * Very useful for array with known dimension. Emulate
 * the original C style array
 * @return The indexed element
 */
dtype & ndarray<dtype>:: operator ()(int co0, ...){
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
//template <typename dtype>
/* *Deprecated*
 * @brief Roll the index in a given axis with periodical boundary condition
 * @param rawind raw index
 * @param axis
 * @return Corresponding value after roll
 */
//dtype ndarray<dtype>::rollval(int rawind, int axis){
//    return head[rollindex(rawind,axis)];
//}
//Other
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
 * @brief Set all the values in the array to a init value
 * @param val   initial value
 */
ndarray<dtype>& ndarray<dtype>:: operator= (dtype val){
    for(int i=0;i<size();i++)
        head[i]=val;
    return *this;
}

#endif //NDARRAY_H
