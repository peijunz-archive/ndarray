#ifndef NDARRAY_H
#define NDARRAY_H
typedef char dtype;
struct ndarray{
    int dim;//dimension of the array
    int flag;
    int *shape;//shape of each dimension
    int *stride;//stride of
    int *shift;//shift of a specific axis
    dtype* head; //head of the array
};
void printarr(ndarray arr);
void setmem(ndarray &arr, int dim, int *shape);
void destroy(ndarray &arr);
//inline int size(ndarray &arr);
void init(ndarray &arr, dtype initval=0);
void init(ndarray &arr, dtype *initarr);
int getval(ndarray &arr, int *index, int flag=0);
int setval(ndarray &arr, int *index, dtype val, int flag=0);
dtype rollval(ndarray arr, int rawind, int axis, bool pn=true);
void place(ndarray arr, double prob=0.5);
#endif //NDARRAY_H
