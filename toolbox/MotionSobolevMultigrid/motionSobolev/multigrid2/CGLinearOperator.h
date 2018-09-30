#ifndef CGLINEAROPERATOR_H
#define CGLINEAROPERATOR_H

template<class T>
class CGLinearOperator {
 public:
  virtual void operator() ( const T *x, T *Ax ) = 0;
};

#endif
