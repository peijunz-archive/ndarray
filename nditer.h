#ifndef nditer_H
#define nditer_H
#include <iostream>
using namespace std;
/**
 * @file nditer.h
 * @author zpj
 * @brief Iterator of sub-index for an ndarray
 * @bug No
 */
/// Iterator of sub-index for an ndarray
template<int D, class size_t=int>
class nditer{
private:
    size_t _shape[D]; ///< Max Bound
public:
    /// An array restoring current sub-index
    size_t _index[D+1];
    nditer(initializer_list<uint> l){
        copy(begin(l), end(l), _shape);
        for(size_t i=0;i<D+1;i++){
            _index[i]=0;
        }
    }
    nditer(const size_t *sh){
        copy(sh, sh+D, _shape);
        for(size_t i=0;i<D+1;i++){
            _index[i]=0;
        }
    }

    /**
     * @brief Go one step further for corresponding raw index
     * @return If the iterator reached the end or not
     */
    bool operator++(){
        _index[D]+=1;
        for(size_t j=D-1;_index[j+1]==_shape[j];j--){
            _index[j+1]=0;
            _index[j]+=1;
        }
        return !_index[0];
	}
    void print(){
        for(int i=0;i<D;i++)
            cout<<_index[i+1]<<", ";
        cout<<endl;
    }
    const size_t * operator()() const{
        return _index+1;
	}
    const size_t * shape() const{
        return _shape;
    }
};

#endif
