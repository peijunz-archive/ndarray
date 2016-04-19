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
 * @todo
 * + broadcasting?
 * + `print()` --> `<<`
 * + Math operations
 */
template <typename dtype>
/// The n-Dimensional C-style array template class
class ndarray{
private:
    int dim;            ///< dimension of the array
    int *shape;         ///< shape of each dimension
    /**
     * @brief stride of every axis.
     *
     * For the purpose of transposition, the naive stride of last index is reserved.
     * So the last stride may be non-zero. Be cautious!
     */
    int *stride;
public:
    dtype *head;        ///< head of the array data
    /// Clear Memory
    void clear(){
        if(dim){
            delete [] shape;
            delete [] stride;
            delete [] head;
            stride=&dim;
        }
    }
    ~ndarray<dtype>(){clear();}

    /// Naive Constructor setting size to 0
    ndarray():dim(0),shape(nullptr),stride(&dim),head(nullptr){}
    /**
     * @brief Initialize an ndarray for a square matrix
     * @param d     dimension
     * @param w     width
     */
    ndarray(int d, int w):
        dim(d),shape(new int[d]),stride(new int[d+1]){
        stride[d]=1;
        for(int i=0;i<d;i++){
            shape[i]=w;
            stride[d-i-1]=stride[d-i]*w;
        }
        head=new dtype[*stride];
    }
    /**
     * @brief A ndarray with a shape given by pointer/array
     * @param d     dimension
     * @param sh    shape array
     *
     * For dynamic construction
     */
    ndarray(int d, int *sh):
        dim(d),shape(new int[d]),stride(new int[d+1]){
        stride[d]=1;
        for(int i=0;i<d;i++){
            shape[i]=sh[i];
            stride[d-i-1]=stride[d-i]*sh[d-i-1];
        }
        head=new dtype[*stride];
    }
    /**
     * @brief Initialize an ndarray with a shape given by an initializer list
     * @param l     initializer list for shape
     *
     * For static construction
     */
    ndarray(initializer_list<int> l):dim(l.size()){
        shape=new int[dim];
        copy(begin(l), end(l), shape);
        stride=new int[dim+1];
        stride[dim]=1;
        for(int i=0;i<dim;i++){
            stride[dim-i-1]=stride[dim-i]*shape[dim-i-1];
        }
        head=new dtype[*stride];
    }

    /**
     * @brief Copy constructor copying the shape only.
     * @param x     The source of ndarray
     * + You can copy the shape at the definition or copy the shape with moving assignment anywhere
     * + If you want to copy the data, use copy assignment!
     *
     * ```cpp
     * ndarray<int> a{3,4};
     * a=1;
     * ndarray<int> b=a, c; // b get the shape of a at definition
     * c=ndarray<int>(a)    // c get the shape of a
     * ```
     */
    ndarray(const ndarray<dtype> & x):
        dim(x.dim),shape(new int[x.dim]),\
        stride(new int[x.dim+1]),\
        head(new dtype[*(x.stride)])
    {
        for(int i=0;i<dim;i++){
            shape[i]=x.shape[i];
            stride[i]=x.stride[i];
        }
        stride[dim]=x.stride[dim];
    }
    /// Move constructor
    ndarray(ndarray<dtype> && x):
        dim(x.dim),shape(x.shape),head(x.head),stride(x.stride){
        x.shape=nullptr;
        x.head=nullptr;
        x.stride=nullptr;
    }
    /**
     * @brief Copy Assignment copying all data
     * @param x     The source of ndarray
     * @todo 形状相同但是权重不同时候的复制，可能需要用到nditer辅助？进行Shape匹配的复制？
     *
     * + 复制空数列等于清空被赋值的数列。
     * + 不能进行size大小不同的复制导致重新分配内存的做法，除非复制到空数列。
     * + 间接方法是先clear再赋值复制。
     * + size相同的时候按照裸指标赋值，无论stride是否相同。
     *
     * ```
     * ndarray<int> c{3,3}, d;
     * c=1
     * d=c;    // Copy the whole data of c
     * ```
     */
    ndarray<dtype> &  operator= (const ndarray<dtype>& x){
        if(*(x.stride)==0){         // Implicitly Clear the array
            clear();
            return *this;
        }
        if(*stride==0){           // After copying, the shape become equal
            *this=ndarray<dtype>(x);
        }
        if(*stride==*(x.stride)){
            for(int i=0;i<*stride;i++){
                head[i]=x.head[i];
            }
        }
        else{
            cerr<<"Unmatched copy! Please clear() first."<<endl;
            exit(0);
        }
        return *this;
    }
    /// Move assignment
    ndarray<dtype> &  operator= (ndarray<dtype>&& x){
        dim=x.dim;
        shape=x.shape;
        head=x.head;
        stride=x.stride;
        x.shape=nullptr;
        x.head=nullptr;
        x.stride=nullptr;
        return *this;
    }
    /**
     * @brief Set all the values in the array to a init value
     * @param val   initial value
     */
    ndarray<dtype>&  operator= (dtype val){
        for(int i=0;i<*stride;i++)
            head[i]=val;
        return *this;
    }

    int size() {return *stride;}      ///< The size of data as an mask of stride[0]
    int _shape(int i){return shape[i];}    ///< Get shape
    int _dim(){return dim;}                ///< Get dimension
    int _stride(int i){return stride[i+1];}///< Get stride

    dtype & operator[](int rawind){return *(head+rawind);}/// Naive indexing
    /**
     * @brief Indexing by an index array
     * @param coo,... the coordinates of an element
     *
     * Especially useful for array with known dimension.
     * @return The indexed element
     */
    dtype & operator()(int *coo){
        int offset=0;
        for(int i=0;i<dim;i++){
            offset+=stride[i+1]*coo[i];
        }
        return head[offset];
    }
    /**
     * @brief Indexing by the argument list of operator()
     * @param co0,... the index of an element
     *
     * Very useful for array with known dimension. Emulate
     * the original C style array
     * @return The indexed element
     */
    dtype & operator ()(int co0, ...){
        va_list l;
        va_start(l, co0);
        int offset=co0*stride[1];
        for(int i=1;i<dim;i++){
            offset+=stride[i+1]*va_arg(l,int);
        }
        va_end(l);
        return head[offset];
    }

    /**
     * @brief Roll the index in a given axis with periodical boundary condition
     * @param rawind raw index
     * @param axis
     * @return Raw index after rolling
     */
    int rollindex(int rawind, int axis){
        int axisind=(rawind%stride[axis])/stride[axis+1];
        if(axis<0){
            axis+=dim;
            if(axisind==0)
                rawind+=stride[axis];
            rawind-=stride[axis+1];
        }
        else{
            if(axisind+1==shape[axis])
                rawind-=stride[axis];
            rawind+=stride[axis+1];
        }
        return rawind;
    }
    /// Transpose between two axis
    void transpose(int i=1, int j=0){
        if(i==j) return;
        int tmp;
        tmp=stride[i+1];
        stride[i+1]=stride[j+1];
        stride[j+1]=tmp;
        tmp=shape[i];
        shape[i]=shape[j];
        shape[j]=tmp;
    }
    /**
     * @brief print the matrix/ndarray
     *
     * + 1D array, printed horizontally
     * + 2D array, printed as an matrix
     * + 3+D array, printed all indices and corresponding values.
     */
    void print(){
        int i,j;
        int ind[dim+1]={0};
        if(dim==0){
            cout<<"NO ARRAY!"<<endl;
        }
        else if (dim==1){
            cout<<"1D array:"<<endl;
            for(i=0;i<*stride;i++)
                cout<<head[i]<<'\t';
            cout<<endl;
            return;
        }
        else if(dim==2){
            for(i=0;i<shape[0];i++){
                for (j=0;j<shape[1];j++){
                    cout<<head[i*stride[1]+j*stride[2]]<<'\t';
                }
                cout<<endl;
            }
            return;
        }
        cout<<"Dimension:"<<dim<<'\t'<<"Size:"<<*stride<<'\n';
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
        for(i=0;i<*stride;i++){
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
};

#endif //NDARRAY_H
