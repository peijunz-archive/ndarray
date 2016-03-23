#ifndef NDARRAY_H
#define NDARRAY_H

typedef char dtype;

struct ndarray{
    int dim;//dimension of the array
    int *shape;//shape of each dimension
    int *stride;//stride[0] is size
    dtype *head; //head of the array
};

void setmem(ndarray &arr, int dim, int *shape);

void destroy(ndarray &arr);

void init(ndarray &arr, dtype initval=0);//Tested

dtype rollval(ndarray arr, int rawind, int axis);

int rollindex(ndarray arr, int rawind, int axis, int &dw);

void place(ndarray arr, double prob=0.5);//Tested

int getval(ndarray &arr, int *index, int flag=0);//Tested

void setval(ndarray &arr, int *index, dtype val, int flag=0);//Tested

void printarr(ndarray &arr);//Tested

#endif //NDARRAY_H
