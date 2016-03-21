#ifndef NDARRAY_H
#define NDARRAY_H
typedef char dtype;
struct ndarray{
    int dim;//dimension of the array
    int flag;
    int size;//dimension of the array
    int *shape;//shape of each dimension
    int *stride;//stride of
    int *shift;//shift of a specific axis
    dtype* head; //head of the array
};
void printarr(ndarray arr);
void setmem(ndarray &array, int dim, int *shape);
void destroy(ndarray &arr);
void init(ndarray &array, dtype initval=0);
void init(ndarray &array, dtype *initarr);
int getval(ndarray &array, int *index, int flag=0);
int setval(ndarray &array, int *index, dtype val, int flag=0);
#endif //NDARRAY_H
