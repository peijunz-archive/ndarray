#ifndef NDARRAY_H
#define NDARRAY_H
#include <stdio.h>

typedef char dtype;

class ndarray{
public:
    int dim;//dimension of the array
    int *shape;//shape of each dimension
    int *stride;//stride of
    dtype *head; //head of the array
    ndarray(int d, int width);
    void reset(int d, int w);
    ndarray(int dim, int *shape);
    ndarray():dim(0),shape(NULL),stride(NULL){}
    ~ndarray();
    inline int size();
    void print();
    int rollindex(int rawind, int axis);
    dtype rollval(int rawind, int axis);
};

void randomize(ndarray &a, double prob);

void homogenize(ndarray &a, dtype initval=0);

#endif //NDARRAY_H
