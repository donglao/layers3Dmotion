#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <math.h>
#include <ostream>

using namespace std;

template <class T> class vector2D {

 public:
   
  T x,y;

 public:

   //Constructors

   vector2D() {;}
   vector2D(T scalar) {x=y=scalar;}
   vector2D(T xx, T yy) {x=xx; y=yy;}

   //Copy constructor which allows casting from vector2D<K> to vector2D<T>

   template<class K> vector2D(const vector2D<K> & V) {x=V.x; y=V.y;}

   //Operations
   void operator = (vector2D V) {x=V.x; y=V.y;}
   void operator = (T scalar) {x=scalar; y=scalar;}

   vector2D operator - () const {return vector2D(-x,-y);}   //Unary -

   vector2D operator + (const vector2D& V) const {return vector2D(x+V.x, y+V.y);}
   vector2D operator - (const vector2D& V) const {return vector2D(x-V.x, y-V.y);}

   void operator += (const vector2D& V) {x+=V.x; y+=V.y;}
   void operator -= (const vector2D& V) {x-=V.x; y-=V.y;}

   T operator * (const vector2D& V) const {return x*V.x + y*V.y;}

   int operator == (vector2D V) {return x==V.x && y==V.y;}
   int operator != (vector2D V) {return x!=V.x || y!=V.y;}

   int operator > (vector2D V)  {return *this * *this > V * V;}
   int operator < (vector2D V)  {return *this * *this < V * V;}
   int operator >= (vector2D V) {return *this * *this >= V * V;}
   int operator <= (vector2D V) {return *this * *this <= V * V;}

   vector2D operator * (T scalar) const {return vector2D(x*scalar,y*scalar);}
   vector2D operator / (T scalar) const {return vector2D(x/scalar,y/scalar);}

   void operator *= (T scalar) {x*=scalar; y*=scalar;}
   void operator /= (T scalar) {x/=scalar; y/=scalar;}

   //Member functions

   double     norm() {return sqrt(x*x+y*y);}
   double   normsq() {return      x*x+y*y ;}
   double   normL1() {return fabs(x)+fabs(y);}
   double normLinf() {return fabs(x)>fabs(y) ? fabs(x) : fabs(y);}


   vector2D rotate(double theta) {return *this.rotateClockwise(theta);}

   vector2D rotateClockwise(double theta)
    {double c=cos(theta),s=sin(theta); return vector2D(x*c+y*s,y*c-x*s);}

   vector2D rotateCounterClockwise(double theta)
    {double c=cos(theta),s=sin(theta); return vector2D(x*c-y*s,y*c+x*s);}


   vector2D rotate90() {return *this.rotateClockwise90();}

   vector2D rotateClockwise90() {return vector2D(y,-x);}

   vector2D rotateCounterClockwise90() {return vector2D(-y,x);}

   //Friend functions

   friend vector2D operator *(T scalar,vector2D V) {return V*scalar;}

   friend double magnitude(vector2D vec) {return vec.norm();}
   friend double      norm(vector2D vec) {return vec.norm();}
   friend double    normsq(vector2D vec) {return vec.normsq();}
   friend double    normL1(vector2D vec) {return vec.normL1();}
   friend double  normLinf(vector2D vec) {return vec.normLinf();}

   friend ostream& operator<< (ostream & os, const vector2D &V)
    {os << '(' << V.x << ',' << V.y << ')'; return os;}
};

#endif
