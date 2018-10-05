#ifndef NDARRAY_H
#define NDARRAY_H
#include <iostream>
#include <algorithm>
#include <initializer_list>
using namespace std;
/**
 * @file ndarray.h
 * @author zpj
 * @brief The n-Dimensional array template for n>1
 *
 * Should have a view class for slicing?
 */

template <typename dtype, int D>
/// The n-Dimensional C-style array template class
class ndarray {
  public:
    using size_t=int32_t;
    using iterator=dtype *;
    using const_iterator=const dtype *;

    ~ndarray() {
        delete [] head;
    }

    ndarray(): head(nullptr) {}

    explicit ndarray(size_t width){
        fill(_shape, _shape+D, width);
        allocate_data();
    }

    explicit ndarray(const size_t *sh){
        copy(sh, sh+D, _shape);
        allocate_data();
    }

    explicit ndarray(initializer_list<size_t> l){
        copy(l.begin(), l.end(), _shape);
        allocate_data();
    }

    template<typename T>
    ndarray(const ndarray<T, D> & x){
        copy(x._shape, x._shape+D, _shape);
        copy(x._stride, x._stride+D+1, _stride);
        head = new dtype[size()];
        copy(x.cbegin(), x.cend(), head);
    }

    ndarray(ndarray<dtype, D> && x):head(x.head){
        copy(x._shape, x._shape+D, _shape);
        copy(x._stride, x._stride+D+1, _stride);
        x.head = nullptr;
    }

    ndarray<dtype, D> & operator= (const ndarray<dtype, D>& x) {
        if(this != &x) {
            ndarray<dtype, D> lhs(x);
            swap(*this, lhs);
        }
        return *this;
    }

    ndarray<dtype, D> & operator= (ndarray<dtype, D>&& x) {
        copy(x._shape, x._shape+D, _shape);
        copy(x._stride, x._stride+D+1, _stride);
        swap(head, x.head);
        return *this;
    }

    ndarray<dtype, D> & operator= (dtype val) {
        fill(begin(), end(), val);
        return *this;
    }


    template <class IndexIterator>
    dtype & operator[](const IndexIterator coo) {
        size_t offset = 0;
        for(int i = 0; i < D-1; i++) {
            offset += coo[i] * _stride[i + 1];
        }
        offset += coo[D-1];
        return head[offset];
    }
    template <class IndexIterator>
    const dtype & operator[](const IndexIterator coo) const{
        size_t offset = 0;
        for(int i = 0; i < D-1; i++) {
            offset += coo[i] * _stride[i + 1];
        }
        offset += coo[D-1];
        return head[offset];
    }

    dtype & operator[](size_t rawind) {
        return head[rawind];
    }
    const dtype & operator[](size_t rawind) const {
        return head[rawind];
    }

    template<int axis>
    size_t adder(size_t last) const {
        return _stride[axis+1] * last;
    }

    template<int axis, typename... Args>
    size_t adder(size_t first, Args... args) const {
        return _stride[axis+1] * first + adder<axis+1>(args...);
    }

    /**
     * @brief Indexing by the argument list of operator()
     * @param args the index of an element
     *
     * + If args are not enough, we fill it with 0
     *
     * @return The indexed element
     */
    template<typename... Args>
    dtype & operator()(Args... args) {
        return head[adder<0>(args...)];
    }
    template<typename... Args>
    const dtype & operator()(Args... args) const {
        return head[adder<0>(args...)];
    }

    ndarray<dtype, D>& operator-() const{
        ndarray<dtype, D> tmp(*this);
        for(size_t i=0;i<size();i++){
            tmp[i]=-head[i];
        }
        return tmp;
    }
    template<int k>
    /**
     * Adding rhs with some constexpr coefficient
     * Also foundation for plus, minus
     */
    ndarray<dtype, D>& plus_k_rhs(const ndarray<dtype, D> &rhs){
        for(size_t i = 0; i < size(); i++)
            head[i] += k*rhs[i];
        return *this;
    }
    ndarray<dtype, D>& operator+=(const ndarray<dtype, D> &rhs){
        return this->plus_k_rhs<1>(rhs);
    }
    ndarray<dtype, D>& operator-=(const ndarray<dtype, D> &rhs){
        return *this->plus_k_rhs<-1>(rhs);
    }
    ///< Not using += and -= directly because of the copying cost
    ndarray<dtype, D>& operator+(const ndarray<dtype, D> &rhs) const {
        ndarray<dtype, D> tmp(*this);
        for(size_t i=0;i<size();i++){
            tmp[i]=head[i]+rhs[i];
        }
        return tmp;
    }
    ndarray<dtype, D>& operator-(const ndarray<dtype, D> &rhs) const {
        ndarray<dtype, D> tmp(*this);
        for(size_t i=0;i<size();i++){
            tmp[i]=head[i]-rhs[i];
        }
        return tmp;
    }
    iterator begin() const {
        return head;
    }
    iterator end() const {
        return head+size();
    }
    const_iterator cbegin() const {
        return head;
    }
    const_iterator cend() const {
        return head+size();
    }
    size_t size() const {return *_stride;}
    int dim() const {return D;}
    size_t shape(int i) const {return _shape[i];}
    const size_t * shape() const {return _shape;}
    size_t stride(int i) const {return _stride[i + 1];}
    const size_t * stride() const {return _stride + 1;}
    const size_t * raw_stride() const {return _stride ;}
  private:
    static_assert(D > 0, "Dimension should be positive integer!");
    size_t __restrict__ _shape[D] = {0};         ///< shape of each dimension
    size_t __restrict__ _stride[D+1] = {0};                  ///< stride of every axis
    dtype* __restrict__ head = nullptr;        ///< head of the array data

    void allocate_data() {
        _stride[D] = 1;
        for(int i = 0; i < D; i++){
            if(_shape[i] <= 0) {
                cerr << "Error: Positive shape needed!" << endl;
                exit(0);
            }
            _stride[D - i - 1] = _stride[D - i] * _shape[D - i - 1];
        }
        head = new dtype[size()];
    }
};

template<bool forward>
/**
* @brief Roll the index in a given axis with periodical boundary condition
* @param rawind raw index
* @param axis
* @param forward   positive/negative
* + For positive axis: `0, 1,..., dim-1`
* + For corresponding negative axis: `-dim, 1-dim,..., -1`
* @return Raw index after rolling
*/
size_t rollindex(size_t rawind, int axis, const size_t* shape, const size_t* raw_stride) {
    size_t axisind = (rawind % raw_stride[axis]) / raw_stride[axis + 1];
    rawind += (2*forward-1)*raw_stride[axis + 1];
    if(axisind == forward*(shape[axis]-1))
        rawind -= (2*forward-1)*raw_stride[axis];
    return rawind;
}
/**
* @brief roll ind for unknow positiveness by using rollindex
*/
template <typename T, int D>
size_t rollind(size_t rawind, int ax, const ndarray<T, D> &arr) {
    if(ax < 0) {
        return rollindex<false>(rawind, ax + D, arr.shape(), arr.raw_stride());
    } else {
        return rollindex<true>(rawind, ax, arr.shape(), arr.raw_stride());
    }
}

template<class dtype>
class matrix: public ndarray<dtype, 2> {
  public:
    using ndarray<dtype, 2>::operator=;
    using ndarray<dtype, 2>::ndarray;
    friend ostream& operator<<(ostream &os, const matrix &m) {
        cout<< "Matrix: "<<m.shape(0)<<"x"<<m.shape(1)<<endl;
        for(int i = 0; i < m.shape(0); i++) {
            for(int j = 0; j < m.shape(1); j++) {
                os << m(i, j) << '\t';
            }
            os << endl;
        }
        return os;
    }
    /**
     * Matrix product (deprecated due to inefficiency)
     */
    matrix<dtype> operator*(matrix<dtype> &rhs) {
        matrix<dtype> c{this->shape(0), rhs.shape(1)};
        c = 0;
        if(this->shape(1) != rhs.shape(0)) {
            cerr << "Dimension not match!" << endl;
            return c;
        }
        for(int i = 0; i < c.shape(0); i++) {
            for(int j = 0; j < c.shape(1); j++) {
                for(int k = 0; k < this->shape(1); k++) {
                    c(i, j) += (*this)(i, k) * rhs(k, j);
                }
            }
        }
        return c;
    }
    int nrow() const {
        return this->shape(0);
    }
    int ncol() const {
        return this->shape(1);
    }
};

#endif //NDARRAY_H
