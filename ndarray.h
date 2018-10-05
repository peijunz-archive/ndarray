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
    /// Constructor an empty ndarray
    ndarray(): head(nullptr) {}
    inline void init_strides(){
        _stride[D] = 1;
        for(int i = 0; i < D; i++)
            _stride[D - i - 1] = _stride[D - i] * _shape[D - i - 1];
    }
    /**
     * @brief Initialize an ndarray for a square matrix
     * @param width
     */
    explicit ndarray(size_t width){
        for(int i = 0; i < D; i++) {
            _shape[i] = width;
        }
        check();
        init_strides();
        head = new dtype[size()];
    }
    /**
     * @brief A ndarray with a shape given by pointer/array
     * @param sh    shape array
     *
     * For dynamic construction
     */
    explicit ndarray(const size_t *sh){
        copy(sh, sh+D, _shape);
        check();
        init_strides();
        head = new dtype[size()];
    }
    /**
     * @brief Initialize an ndarray with a shape given by an initializer list
     * @param l     initializer list for shape
     *
     * For static construction
     */
    explicit ndarray(initializer_list<size_t> l){
        copy(l.begin(), l.end(), _shape);
        check();
        init_strides();
        head = new dtype[size()];
    }

    /**
     * @brief Copy constructor copying the shape only.
     * @param x     The source of ndarray
     */
    template<typename T>
    ndarray(const ndarray<T, D> & x){
        copy(x._shape, x._shape+D, _shape);
        copy(x._stride, x._stride+D+1, _stride);
        head = new dtype[size()];
        copy(x.cbegin(), x.cend(), head);
    }
    /// Move constructor
    ndarray(ndarray<dtype, D> && x):head(x.head){
        copy(x._shape, x._shape+D, _shape);
        copy(x._stride, x._stride+D+1, _stride);
        x.head = nullptr;
    }
    /**
     * @brief Copy Assignment copying all data
     * @param x     The source of ndarray
     */
    ndarray<dtype, D> & operator= (const ndarray<dtype, D>& x) {
        if(this != &x) {
            ndarray<dtype, D> lhs(x);
            swap(*this, lhs);
        }
        return *this;
    }
    /// Move assignment
    ndarray<dtype, D> & operator= (ndarray<dtype, D>&& x) {
        copy(x._shape, x._shape+D, _shape);
        copy(x._stride, x._stride+D+1, _stride);
        swap(head, x.head);
        return *this;
    }
    /**
     * @brief Set all the values in the array to a init value
     * @param val   initial value
     */
    ndarray<dtype, D> & operator= (dtype val) {
        for(size_t i = 0; i < size(); i++)
            head[i] = val;
        return *this;
    }
    size_t size() const {
        return *_stride;
    }
    int dim() const {
        return D;
    }
    size_t shape(int i) const {
        return _shape[i];
    }
    const size_t * shape() const {
        return _shape;
    }
    size_t stride(int i) const {
        return _stride[i + 1];
    }
    const size_t * stride() const {
        return _stride + 1;
    }

    /**
     * @brief Indexing by an index array
     * @param coo   the coordinates of an element
     * @return The indexed element
     */
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
    /**
     * @brief Naive indexing using raw index
     * @param rawind    raw index
     *
     * + Good efficiency!
     * + If dim==1, use [] rather than ()
     */
    dtype & operator[](size_t rawind) {
        return *(head + rawind);
    }
    const dtype & operator[](size_t rawind) const {
        return *(head + rawind);
    }


    size_t adder(size_t k, size_t last) const {
        return last;//k=1
    }
    template<typename... Args>
    size_t adder(size_t k, size_t first, Args... args) const {
        return _stride[k] * first + adder(k + 1, args...);
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
    dtype & operator()(Args... args) {
        return *(head + adder(1, args...));
    }
    template<typename... Args>
    const dtype & operator()(size_t arg0, Args... args) const {
        return *(head + adder(1, args...));
    }

    template<bool dir>
    /**
     * @brief Roll the index in a given axis with periodical boundary condition
     * @param rawind raw index
     * @param axis
     * @param dir   positive/negative
     * + For positive axis: `0, 1,..., dim-1`
     * + For corresponding negative axis: `-dim, 1-dim,..., -1`
     * @return Raw index after rolling
     */
    size_t rollindex(size_t rawind, int axis) const {
        size_t axisind = (rawind % _stride[axis]) / _stride[axis + 1];
        rawind += (2*dir-1)*_stride[axis + 1];
        if(axisind == dir*(_shape[axis]-1))
            rawind -= (2*dir-1)*_stride[axis];
        return rawind;
    }
    /**
     * @brief roll ind for unknow positiveness by using rollindex
     */
    size_t rollind(size_t rawind, int ax) const {
        if(ax < 0) {
            return rollindex<false>(rawind, ax + D);
        } else {
            return rollindex<true>(rawind, ax);
        }
    }

    /// Transpose between two axis
    void transpose(int i = 1, int j = 0) {
        swap(_stride[i + 1], _stride[j + 1]);
        swap(_shape[i + 1], _shape[j + 1]);
    }
    int size_attached() const { ///< empty condition
        int s = sizeof(ndarray<dtype, D>);
        s += sizeof(size_t) * (2 * D + 1);
        return s;
    }
    ndarray<dtype, D>& operator-() const{
        ndarray<dtype, D> tmp(*this);
        for(size_t i=0;i<size();i++){
            tmp[i]=-head[i];
        }
        return *this;
    }
    template<int N>
    /**
     * Adding rhs with some constexpr coefficient
     * Also foundation for plus, minus
     */
    ndarray<dtype, D>& pm(const ndarray<dtype, D> &rhs){
        for(size_t i = 0; i < size(); i++)
            head[i] += N*rhs[i];
        return *this;
    }
    ndarray<dtype, D>& operator+=(const ndarray<dtype, D> &rhs){
        return this->pm<1>(rhs);
    }
    ndarray<dtype, D>& operator-=(const ndarray<dtype, D> &rhs){
        return *this->pm<-1>(rhs);
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
  private:
    static_assert(D > 0, "Dimension should be positive!");
    size_t __restrict__ _shape[D] = {0};         ///< shape of each dimension
    /**
     * @brief stride of every axis.
     *
     * For the purpose of transposition,
     * the naive stride of last index is reserved.
     * So the last stride may be non-zero. Be cautious!
     */
    size_t _stride[D+1] = {0};
    dtype* __restrict__ head = nullptr;        ///< head of the array data
    void check() {
        for(int i = 0; i < D; i++) {
            if(_shape[i] <= 0) {
                cerr << "Error: Positive shape needed!" << endl;
                exit(0);
            }
        }
    }
};

template<class dtype>
/**
 * Two dimensional array---matrix
 */
class matrix: public ndarray<dtype, 2> {
  public:
    using ndarray<dtype, 2>::operator=;
    using ndarray<dtype, 2>::ndarray;
    void print() {
        for(int i = 0; i < this->shape(0); i++) {
            for(int j = 0; j < this->shape(1); j++) {
                cout << (*this)(i, j) << '\t';
            }
            cout << endl;
        }
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
