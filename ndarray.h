#ifndef NDARRAY_H
#define NDARRAY_H
#include <iostream>
#include <vector>
#include <initializer_list>
#define CPP11SUPPORT
using namespace std;
/**
 * @file ndarray.h
 * @author zpj
 * @brief The n-Dimensional array template for n>1
 * @todo
 * + broadcasting?
 * + Math operations
 * + Make dimension a template argument
 *
 * Should have a view class for slicing?
 */
template <typename dtype, uint D>
/// The n-Dimensional C-style array template class
class ndarray {
  private:
    uint *_shape;         ///< shape of each dimension
    /**
     * @brief stride of every axis.
     *
     * For the purpose of transposition,
     * the naive stride of last index is reserved.
     * So the last stride may be non-zero. Be cautious!
     */
    uint *_stride;
    void check() {
        for(uint i = 0; i < D; i++) {
            if(_shape[i] <= 0) {
                cerr << "Error: Positive shape needed!" << endl;
                exit(0);
            }
        }
    }

  public:
    dtype *head;        ///< head of the array data
    ~ndarray() {
        delete [] head;
        delete [] _shape;
        delete [] _stride;
    }
    /// Naive Constructor setting size to 0
    ndarray(): _shape(new uint[D]), _stride(new uint[D + 1]), head(nullptr) {
        *_stride = 0;
    }
    /**
     * @brief Initialize an ndarray for a square matrix
     * @param width
     */
    ndarray(uint width): _shape(new uint[D]), _stride(new uint[D + 1]) {
        _stride[D] = 1;
        for(uint i = 0; i < D; i++) {
            _shape[i] = width;
            _stride[D - i - 1] = _stride[D - i] * width;
        }
        check();
        head = new dtype[*_stride];
    }
    /**
     * @brief A ndarray with a shape given by pointer/array
     * @param sh    shape array
     *
     * For dynamic construction
     */
    ndarray(uint *sh): _shape(new uint[D]), _stride(new uint[D + 1]) {
        _stride[D] = 1;
        for(uint i = 0; i < D; i++) {
            _shape[i] = sh[i];
            _stride[D - i - 1] = _stride[D - i] * sh[D - i - 1];
        }
        check();
        head = new dtype[*_stride];
    }
#ifdef CPP11SUPPORT
    /**
     * @brief Initialize an ndarray with a shape given by an initializer list
     * @param l     initializer list for shape
     *
     * For static construction
     */
    ndarray(initializer_list<uint> l): _shape(new uint[D]), _stride(new uint[D + 1]) {
        if(l.size() != D)
            exit(0);
        copy(begin(l), end(l), _shape);
        _stride[D] = 1;
        for(uint i = 0; i < D; i++) {
            _stride[D - i - 1] = _stride[D - i] * _shape[D - i - 1];
        }
        check();
        head = new dtype[*_stride];
    }
#endif

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
    ndarray(const ndarray<T, D> & x): _shape(new uint[D]), _stride(new uint[D + 1]) {
        for(uint i = 0; i < D; i++) {
            _shape[i] = x.shape(i);
            _stride[i + 1] = x.stride(i);
        }
        _stride[0] = x.size();
        head = new dtype[x.size()];
    }
    /// Move constructor
    ndarray(ndarray<dtype, D> && x): _shape(x._shape), head(x.head), _stride(x._stride) {
        x._shape = nullptr;
        x.head = nullptr;
        x._stride = nullptr;
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
    ndarray<dtype, D> & operator= (const ndarray<dtype, D>& x) {
        if(this != &x) {
            if(*(x._stride) == 0) {
                delete [] head;
                return *this;
            }
            if(*_stride == 0) {
                *this = ndarray<dtype, D>(x);
            }
            if(*_stride == *(x._stride)) {
                for(uint i = 0; i < *_stride; i++) {
                    head[i] = x.head[i];
                }
            } else {
                cerr << "Unmatched copy! Please clear() first." << endl;
                exit(0);
            }
        }
        return *this;
    }
    /// Move assignment
    ndarray<dtype, D> & operator= (ndarray<dtype, D>&& x) {
        _shape = x._shape;
        head = x.head;
        _stride = x._stride;
        x._shape = nullptr;
        x.head = nullptr;
        x._stride = nullptr;
        return *this;
    }
    /**
     * @brief Set all the values in the array to a init value
     * @param val   initial value
     */
    ndarray<dtype, D> & operator= (dtype val) {
        for(uint i = 0; i < *_stride; i++)
            head[i] = val;
        return *this;
    }

    uint size() const {
        return *_stride;
    }
    uint dim() const {
        return D;
    }
    uint shape(uint i) const {
        return _shape[i];
    }
    const uint * shape() const {
        return _shape;
    }
    uint stride(uint i) const {
        return _stride[i + 1];
    }
    const uint * stride() const {
        return _stride + 1;
    }

    /**
     * @brief Naive indexing using raw index
     * @param rawind    raw index
     *
     * + Good efficiency!
     * + If dim==1, use [] rather than ()
     */
    dtype & operator[](uint rawind) {
        return *(head + rawind);
    }

    /**
     * @brief Indexing by an index array
     * @param coo   the coordinates of an element
     * @return The indexed element
     */
    dtype & operator()(uint *coo) {
        uint offset = 0;
        for(uint i = 0; i < D; i++) {
            offset += coo[i] * _stride[i + 1];
        }
        return head[offset];
    }

    uint adder(uint *p, uint last) {
        return *p * last;
    }
    template<typename... Args>
    uint adder(uint *p, uint first, Args... args) {
        return *p * first + adder(p + 1, args...);
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
        return *(head + adder(_stride + 1, args...));
    }

    /**
     * @brief Roll the index in a given axis with periodical boundary condition
     * @param rawind raw index
     * @param axis
     * @param dir   positive/negative
     * + For positive axis: `0, 1,..., dim-1`
     * + For corresponding negative axis: `-dim, 1-dim,..., -1`
     * @return Raw index after rolling
     */
    uint rollindex(uint rawind, int axis, int dir = 1) {
        uint axisind = (rawind % _stride[axis]) / _stride[axis + 1];
        if(dir == 1) {
            rawind += _stride[axis + 1];
            if(axisind + 1 == _shape[axis])
                rawind -= _stride[axis];
        } else {
            rawind -= _stride[axis + 1];
            if(axisind == 0)
                rawind += _stride[axis];
        }
        return rawind;
    }
    /**
     * @brief roll ind for unknow positiveness by using rollindex
     */
    uint rollind(int rawind, int ax) {
        if(ax < 0) {
            return rollindex(rawind, ax + D, -1);
        } else {
            return rollindex(rawind, ax, 1);
        }
    }

    /// Transpose between two axis
    void transpose(uint i = 1, uint j = 0) {
        if(i == j) return;
        uint tmp;
        tmp = _stride[i + 1];
        _stride[i + 1] = _stride[j + 1];
        _stride[j + 1] = tmp;
        tmp = _shape[i];
        _shape[i] = _shape[j];
        _shape[j] = tmp;
    }
    uint size_attached() {
        uint s = sizeof(ndarray<dtype, D>);
        s += sizeof(uint) * (2 * D + 1);
        return s;
    }
    vector<uint> subind(uint ind) {
        vector<uint> t;
        t.reserve(D);
        for(uint i = D - 1; i >= 0; i--) {
            t[i] = ind % _stride[i];
            ind /= _stride[i];
        }
        return t;
    }
};
#endif //NDARRAY_H
