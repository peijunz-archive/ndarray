#ifndef NDARRAY_H
#define NDARRAY_H
#include <iostream>
#include <initializer_list>
using namespace std;
/**
 * @file ndarray.h
 * @author zpj
 * @brief The n-Dimensional array template for n>1
 * @todo
 * + broadcasting?
 * + Math operations
 *
 * Should have a view class for slicing?
 */
template <typename dtype>
/// The n-Dimensional C-style array template class
class ndarray{
private:
    int _dim;            ///< dimension of the array
    int *_shape;         ///< shape of each dimension
    /**
     * @brief stride of every axis.
     *
     * For the purpose of transposition,
     * the naive stride of last index is reserved.
     * So the last stride may be non-zero. Be cautious!
     */
    int *_stride;
    void check(){
        if(_dim<=0){
            cerr<<"Error: Positive dimension needed!"<<endl;
            exit(0);
        }
        for(int i=0;i<_dim;i++){
            if(_shape[i]<=0){
                cerr<<"Error: Positive shape needed!"<<endl;
                exit(0);
            }
        }
    }

public:
    dtype *head;        ///< head of the array data
    /// Clear Memory
    void clear(){
        if(_dim){
            delete [] _shape;
            delete [] _stride;
            delete [] head;
            _stride=&_dim;
        }
    }
    ~ndarray(){clear();}

    /// Naive Constructor setting size to 0
    ndarray():_dim(0),_shape(nullptr),_stride(&_dim),head(nullptr){}
    /**
     * @brief Initialize an ndarray for a square matrix
     * @param d     dimension
     * @param w     width
     */
    ndarray(int d, int w):
        _dim(d),_shape(new int[d]),_stride(new int[d+1]){
        _stride[d]=1;
        for(int i=0;i<d;i++){
            _shape[i]=w;
            _stride[d-i-1]=_stride[d-i]*w;
        }
        check();
        head=new dtype[*_stride];
    }
    /**
     * @brief A ndarray with a shape given by pointer/array
     * @param d     dimension
     * @param sh    shape array
     *
     * For dynamic construction
     */
    ndarray(int d, int *sh):
        _dim(d),_shape(new int[d]),_stride(new int[d+1]){
        _stride[d]=1;
        for(int i=0;i<d;i++){
            _shape[i]=sh[i];
            _stride[d-i-1]=_stride[d-i]*sh[d-i-1];
        }
        check();
        head=new dtype[*_stride];
    }
    /**
     * @brief Initialize an ndarray with a shape given by an initializer list
     * @param l     initializer list for shape
     *
     * For static construction
     */
    ndarray(initializer_list<int> l):_dim(l.size()){
        _shape=new int[_dim];
        copy(begin(l), end(l), _shape);
        _stride=new int[_dim+1];
        _stride[_dim]=1;
        for(int i=0;i<_dim;i++){
            _stride[_dim-i-1]=_stride[_dim-i]*_shape[_dim-i-1];
        }
        check();
        head=new dtype[*_stride];
    }

    /**
     * @brief Copy constructor copying the shape only.
     * @param x     The source of ndarray
     * + You can copy the shape at the definition or copy the
     *   shape with moving assignment anywhere
     * + If you want to copy the data, use copy assignment!
     *
     * ```cpp
     * ndarray<int> a{3,4};
     * a=1;
     * ndarray<int> b=a, c; // b get the shape of a at definition
     * c=ndarray<int>(a)    // c get the shape of a
     * ```
     */
    template<typename T>
    ndarray(const ndarray<T> & x):
        _dim(x.dim()),_shape(new int[x.dim()]),\
        _stride(new int[x.dim()+1]),\
        head(new dtype[x.size()])
    {
        for(int i=0;i<_dim;i++){
            _shape[i]=x.shape(i);
            _stride[i+1]=x.stride(i);
        }
        _stride[0]=x.size();
    }
    /// Move constructor
    ndarray(ndarray<dtype> && x):
        _dim(x._dim),_shape(x._shape),head(x.head),_stride(x._stride){
        x._shape=nullptr;
        x.head=nullptr;
        x._stride=nullptr;
    }
    /**
     * @brief Copy Assignment copying all data
     * @param x     The source of ndarray
     * @todo If the shape is the same while stride is not,
     *   how to copy? Using nditer?
     *
     * + Copy an empty x array will clear() the array
     * + Copy to an empty array will deep copy the array x
     * + Except for above conditions, copy between arrays
     *   of different size is an error
     * + Indirective method is clear first and then copy
     * + If size is equal, copy according to raw index,
     *   regardless of the shape, dim, stride
     *
     * ```
     * ndarray<int> c{3,3}, d;
     * c=1
     * d=c;    // Copy the whole data of c
     * ```
     */
    ndarray<dtype> & operator= (const ndarray<dtype>& x){
        if(this != &x){
            if(*(x._stride)==0){
                clear();
                return *this;
            }
            if(*_stride==0){
                *this=ndarray<dtype>(x);
            }
            if(*_stride==*(x._stride)){
                for(int i=0;i<*_stride;i++){
                    head[i]=x.head[i];
                }
            }
            else{
                cerr<<"Unmatched copy! Please clear() first."<<endl;
                exit(0);
            }
        }
        return *this;
    }
    /// Move assignment
    ndarray<dtype> & operator= (ndarray<dtype>&& x){
        _dim=x._dim;
        _shape=x._shape;
        head=x.head;
        _stride=x._stride;
        x._shape=nullptr;
        x.head=nullptr;
        x._stride=nullptr;
        return *this;
    }
    /**
     * @brief Set all the values in the array to a init value
     * @param val   initial value
     */
    ndarray<dtype> & operator= (dtype val){
        for(int i=0;i<*_stride;i++)
            head[i]=val;
        return *this;
    }

    int size() const {return *_stride;}
    int dim() const {return _dim;}
    int shape(int i) const {return _shape[i];}
    const int * shape() const {return _shape;}
    int stride(int i) const{return _stride[i+1];}
    const int * stride() const{return _stride+1;}

    /**
     * @brief Naive indexing using raw index
     * @param rawind    raw index
     *
     * + Good efficiency!
     * + If dim==1, use [] rather than ()
     */
    dtype & operator[] (int rawind){return *(head+rawind);}

    template<int N>
    /**
     * @brief dot product for small N
     * @param coo   the coordinates of an element
     * @return
     */
    int dot(int *coo){
        int offset=0;
        for(int i=0;i<N;i++){
            offset+=coo[i]*_stride[i+1];
        }
        return offset;
    }
    /**
     * @brief Indexing by an index array
     * @param coo   the coordinates of an element
     *
     * + Especially useful for array with unknown dimension.
     * + The speed is about 0.6+ that of operator ()(Args... args)
     * + If _dim==1, use [] directly rather than ()
     *
     * @return The indexed element
     */
    dtype & operator()(int *coo){
        if(_dim==2) return head[dot<2>(coo)];
        if(_dim==3) return head[dot<3>(coo)];
        if(_dim==4) return head[dot<4>(coo)];
        if(_dim==5) return head[dot<5>(coo)];
        int offset=0;
        for(int i=0;i<_dim;i++){
            offset+=coo[i]*_stride[i+1];
        }
        return head[offset];
    }

    int adder(int *p,int last){
      return *p*last;
    }

    template<typename... Args>
    int adder(int *p, int first, Args... args) {
      return *p*first+adder(p+1, args...);
    }

    /**
     * @brief Indexing by the argument list of operator()
     * @param args the index of an element
     *
     * + Very useful for array with known dimension.
     * + Emulate the original C style array
     * + If _dim==1, use [] directly rather than ()
     * + If args are not enough, suppose we fill it with 0
     *
     * @return The indexed element
     */
    template<typename... Args>
    dtype & operator()(Args... args){
        return *(head+adder(_stride+1, args...));
    }

    /**
     * @brief Roll the index in a given axis with periodical boundary condition
     * @param rawind raw index
     * @param axis
     * + For positive axis: `0, 1,..., dim-1`
     * + For corresponding negative axis: `-dim, 1-dim,..., -1`
     * @return Raw index after rolling
     */
    int rollindex(int rawind, int axis, int dir=1){
        int axisind=(rawind%_stride[axis])/_stride[axis+1];
//        cout<<"rollindex"<<' '<<rawind;
        if(dir==1){
            rawind+=_stride[axis+1];
            if(axisind+1==_shape[axis])
                rawind-=_stride[axis];
        }
        else{
            rawind-=_stride[axis+1];
            if(axisind==0)
                rawind+=_stride[axis];
        }
//        cout<<"->"<<rawind<<" ax "<<axis<<" dir "<<dir<<endl;
        return rawind;
    }
    int rollind(int rawind, int ax, int D){
        if(ax<0){
            return rollindex(rawind, ax+D, -1);
        }
        else{
            return rollindex(rawind, ax, 1);
        }
    }

    /// Transpose between two axis
    void transpose(int i=1, int j=0){
        if(i==j) return;
        int tmp;
        tmp=_stride[i+1];
        _stride[i+1]=_stride[j+1];
        _stride[j+1]=tmp;
        tmp=_shape[i];
        _shape[i]=_shape[j];
        _shape[j]=tmp;
    }
    int size_attached(){
        int s=sizeof(ndarray<dtype>);
        if(_dim){
            s+=sizeof(int)*(2*_dim+1);
        }
        return s;
    }
};
template<typename T>
/**
 * @brief Print a vector for the last dimension of an ndarray to stream
 * @param out   stream
 * @param x     ndarray
 * @param start
 */
void printvec(ostream &out, ndarray<T> &x, const T *start){
    for(int i=0;i<x.shape(x.dim()-1);i++)
        out<<*(start+i*x.stride(x.dim()-1))<<'\t';
    out<<"\n";
}
template<typename T>
/**
 * @brief Print a matrix for the last two dimensions of an ndarray to stream
 * @param out   stream
 * @param x     ndarray
 * @param start
 */
void printmat(ostream &out, ndarray<T> &x, const T *start){
    out<<"\t ";
    for(int i=0;i<x.shape(x.dim()-1);i++) out<<i<<'\t';
    out<<"\n\t|----------------------------------\n";
    for(int i=0;i<x.shape(x.dim()-2);i++){
        out<<i<<"\t|";
        printvec(out, x, start+i*x.stride(x.dim()-2));
    }
    out<<'\n';
}

template<typename T>
/**
 * @brief Print the ndarray
 *
 * + 1 D array, printed horizontally
 * + 2 D array, printed as an matrix
 * + 3+D array, printed as Many Matrics.
 */
void printsubmat(ostream &out, ndarray<T> &x){
    if(x.dim()==0){
        cerr<<"Error: print an empty array"<<endl;
        return;
    }
    else if(x.dim()==1){
        out<<"Size:"<<x.size()<<"\t|";
        printvec(out, x, x.head);
    }
    else{
        out<<"Dimension:"<<x.dim()<<'\t';
        out<<"\nShape\t ";
        for(int i=0;i<x.dim();i++)
            out<<x.shape(i)<<"\t";
        out<<"\n";
        if(x.dim()==2){
            printmat(out, x, x.head);
        }
        else{
            int ind[x.dim()-1]={0};
            int step=x.stride(x.dim()-3), num=x.size()/step;
            for(int i=0;i<num;i++){
                for(int j=x.dim()-3;ind[j+1]==x.shape(j);j--){
                    ind[j+1]=0;
                    ind[j]+=1;
                }
                out<<"Matrix No. "<<i<<": (";
                for(int j=1;j<x.dim()-1;j++)
                    out<<ind[j]<<", ";
                out<<"\b\b) \n";
                printmat(out, x, x.head+i*step);
                ind[x.dim()-2]+=1;
            }
        }
    }
}
template <typename T>
ostream & operator<<(ostream &out, ndarray<int> &x){
    int ind[x.dim()+1]={0};
    for(int i=0;i<x.size();i++){
        for(int j=x.dim()-1;ind[j+1]==x.shape(j);j--){
            ind[j+1]=0;
            ind[j]+=1;
        }
        out<<"Matrix No. "<<i<<": (";
        for(int j=1;j<x.dim()+1;j++)
            out<<ind[j]<<", ";
        out<<"\b\b) :";
        out<<x[i];
        out<<'\n';
        ind[x.dim()]+=1;
    }
    return out;
}

ostream & operator<<(ostream &out, ndarray<int> &x){
    printsubmat(out, x);
    return out;
}
ostream & operator<<(ostream &out, ndarray<double> &x){
    printsubmat(out, x);
    return out;
}
ostream & operator<<(ostream &out, ndarray<char> &x){
    printsubmat(out, x);
    return out;
}
#endif //NDARRAY_H
