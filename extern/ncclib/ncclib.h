////////////////////////////////////////////////////////////////////////////////
//NCCLIB	Naor's C collection library
//
// Contents:
// 	0) Platform specific macros and #includes, e.g. windows.h and Microsoft
// 	compiler macros from the CRT.
//
// 	1) Modified nr3.h. This is a modified version of the Numerical Recipes
// 	   third edition main dependency. To simplify transportability, I decided
// 	   to copy nr3.h and other useful nr routines into ncclib.cpp, rather than
// 	   #include their respective .h files. Since nr doesn't bother with forward
// 	   declarations I had to create those for all nr3 routines and put them
// 	   in here. Note that templated items also have their definition here, and
// 	   not in the corresponding .cpp.
// 	   The contents of nr3.h are surrounded in a nr3 namespace.
//
// 	  1a) Forward declarations for the original nr3.h, inside an nr3 namespace.
// 	      These include a few useful macro names for inline functions (MAX,MIN,
// 	      etc.), an error macro (throw), a quiet NaN variable, and an efficient
// 	      implementation of templated vector and matrix classes.
// 	  1b) Froward declarations for the contents of ran.h, also inside the nr3
// 	      namespace. Much better than the built-in rand(). The usage is:
// 		nr3::Ran ran(seed);
// 		ran.int64(); for a random uint64
// 		ran.doub();  for a random double in (0,1)
// 		1 + ran.int64() % (n-1) for a random uint64 between 1 and n
// 	  1c) Forward declarations for the contents of gamma.h.
// 	  1d) Forward declarations for the contents of deviates.h. Some useful
// 	      probability distributions.
// 	      Usage:
// 		nr3::Normaldev norm(mu.sig,seed);
// 		norm.dev(); for a normal deviate (copied from deviates.h)
//
// 	2) WHOMAI(). An expanded "hello world" message with environment info.
//
// 	3) NCC__ERROR() and NCC__WARNING(). Simple macro printf and exit(1);
//
// 	4) Some file I/O functions.
// 	  4a) GetStrPropertyFromINIFile() and GetIntPropertyFromINIFile()
// 	      Port of GetPrivateProfileString() and GetPrivateProfileInt() from the
// 	      Windows API.
// 	  4b) QDReadMatrixFromFile(). A Quick and Dirty matrix reader using c-style
// 	      arrays and fscans. Doesn't work with gcc. Only useful for demonstrating
// 	      c-style file I/O.
// 	  4c) RSReadMatrixFromFile(). A reasonably safe matrix reader using c++ STL.
// 	      But it doesn't work with gcc. Better to use Load instead anyway.
//    4d) load(). A reasonably safe matrix reader using nr3 matrix types. This
//        is the recommended function for the basic task of reading in data
//        in a rectangular array.
//    4e) logEntry(). A pedestrian convenience function that appends one line to an
//        ascii file.
// 	
//  	5) VEC3F and VEC3D. Rudimentary public access Cartesian 3-vector class.
//
//	6) RGB. A namespace for names colors based on the X11 color spec.
// 	
// Author: Me (Naor Movshovitz)
////////////////////////////////////////////////////////////////////////////////

#ifndef _NCCLIB_
#define _NCCLIB_

// System specific definitions
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#ifdef _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES    //
#undef _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES    //
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1 // use CRT secure versions
#endif
#endif

#ifdef __linux__
#include <unistd.h>
#define _getcwd getcwd
#endif

#ifdef __APPLE__
#include <unistd.h>
#define _getcwd getcwd
#endif

// All the common #includes
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <cstring>
#include <ctype.h>
using namespace std;

// All nr3 entities are part of the nr3 namespace
namespace nr3 {

// Macro-like inline functions (definition in header!)
template<class T> inline T SQR(const T a) {return a*a;}
template<class T> inline const T &MAX(const T &a, const T &b)	{return b > a ? (b) : (a);}
inline float MAX(const double &a, const float &b) 		{return b > a ? (b) : float(a);}
inline float MAX(const float &a, const double &b) 		{return b > a ? float(b) : (a);}
template<class T> inline const T &MIN(const T &a,const T &b)	{return b < a ? (b) : (a);}
inline float MIN(const double &a, const float &b)		{return b < a ? (b) : float(a);}
inline float MIN(const float &a, const double &b)		{return b < a ? float(b) : (a);}
template<class T> inline T SIGN(const T &a, const T &b)		{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
inline float SIGN(const float &a, const double &b)		{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
inline float SIGN(const double &a, const float &b)		{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}
template<class T> inline void SWAP(T &a, T &b)			{T dum=a; a=b; b=dum;}

// Vector and Matrix Classes
// NRvector definitions
template <class T> class NRvector {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRvector();
	explicit NRvector(int n);	// Zero-based array
	NRvector(int n, const T &a);	// Initialize to constant value
	NRvector(int n, const T *a);	// Initialize to array
	NRvector(const NRvector &rhs);	// Copy constructor
	NRvector & operator=(const NRvector &rhs);	// Assignment
	typedef T value_type; // Make T available externally
	inline T & operator[](const int i);	// i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	~NRvector();
};

template <class T> NRvector<T>::NRvector() : nn(0), v(NULL) {}
template <class T> NRvector<T>::NRvector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}
template <class T> NRvector<T>::NRvector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = a;
}
template <class T> NRvector<T>::NRvector(int n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = *a++;
}
template <class T> NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
	for(int i=0; i<nn; i++) v[i] = rhs[i];
}
template <class T> NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != NULL) delete [] (v);
			nn=rhs.nn;
			v= nn>0 ? new T[nn] : NULL;
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}
template <class T> inline T & NRvector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}
template <class T> inline const T & NRvector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}
template <class T> inline int NRvector<T>::size() const
{
	return nn;
}
template <class T> void NRvector<T>::resize(int newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}
template <class T> void NRvector<T>::assign(int newn, const T& a)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i=0;i<nn;i++) v[i] = a;
}
template <class T> NRvector<T>::~NRvector()
{
	if (v != NULL) delete[] (v);
}
// end of NRvector definitions
// NRmatrix definitions
template <class T> class NRmatrix {
private:
	int nn;
	int mm;
	T **v;
public:
	NRmatrix();
	NRmatrix(int n, int m);			// Zero-based array
	NRmatrix(int n, int m, const T &a);	// Initialize to constant
	NRmatrix(int n, int m, const T *a);	// Initialize to array
	NRmatrix(const NRmatrix &rhs);		// Copy constructor
	NRmatrix & operator=(const NRmatrix &rhs);	// Assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	// Subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	~NRmatrix();
};
template <class T> NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}
template <class T> NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}
template <class T> NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}
template <class T> NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}
template <class T> NRmatrix<T>::NRmatrix(const NRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	int i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}
template <class T> NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j,nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}
template <class T> inline T* NRmatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T> inline const T* NRmatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}
template <class T> inline int NRmatrix<T>::nrows() const
{
	return nn;
}
template <class T> inline int NRmatrix<T>::ncols() const
{
	return mm;
}
template <class T> void NRmatrix<T>::resize(int newn, int newm)
{
	int i,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
}
template <class T> void NRmatrix<T>::assign(int newn, int newm, const T& a)
{
	int i,j,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}
template <class T> NRmatrix<T>::~NRmatrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}
// end of NRmatrix definitions
// NRMat3d definitions
template <class T> class NRMat3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	NRMat3d();
	NRMat3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~NRMat3d();
};
template <class T> NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(NULL) {}
template <class T> NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
}
template <class T> inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}
template <class T> inline const T* const * NRMat3d<T>::operator[](const int i) const
{
	return v[i];
}
template <class T> inline int NRMat3d<T>::dim1() const
{
	return nn;
}
template <class T> inline int NRMat3d<T>::dim2() const
{
	return mm;
}
template <class T> inline int NRMat3d<T>::dim3() const
{
	return kk;
}
template <class T> NRMat3d<T>::~NRMat3d()
{
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}
// end of NRMat3d definitions

// basic type names (redefine if your bit lengths don't match)
typedef int Int; // 32 bit integer
typedef unsigned int Uint;
#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif
typedef char Char; // 8 bit integer
typedef unsigned char Uchar;
typedef double Doub; // default floating type
typedef long double Ldoub;
typedef complex<double> Complex; // default complex type
typedef bool Bool;

// vector types
typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;
typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;
typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;
typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;
typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;
typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;
typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;
typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;
typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;
typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;
typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types
typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;
typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;
typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;
typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;
typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;
typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;
typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;
typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types
typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;

// End of original nr3.h
// Begin contents of ran.h

struct Ran {
	Ullong u,v,w;
	Ran(Ullong j);
	Ullong int64();
	Doub doub();
	Uint int32();
};

struct Ranq1 {
	Ullong v;
	Ranq1(Ullong j);
	Ullong int64();
	Doub doub();
	Uint int32();
};

struct Ranq2 {
	Ullong v,w;
	Ranq2(Ullong j);
	Ullong int64();
	Doub doub();
	Uint int32();
};

struct Ranhash {
	Ullong int64(Ullong u);
	Uint int32(Ullong u);
	Doub doub(Ullong u);
};

struct Ranbyte {
	Int s[256],i,j,ss;
	Uint v;
	Ranbyte(Int u);
	unsigned char int8();
	Uint int32();
	Doub doub();
};

struct Ranfib {
	Doub dtab[55], dd;
	Int inext, inextp;
	Ranfib(Ullong j);
	Doub doub();
	unsigned long int32();
};

struct Ranlim32 {
	Uint u,v,w1,w2;
	Ranlim32(Uint j);
	Uint int32();
	Doub doub();
	Doub truedoub();
};

// End of Ran.h.
// Begin contents of gamma.h
Doub gammln(const Doub xx);
Doub factrl(const Int n);
Doub factln(const Int n);
Doub bico(const Int n, const Int k);
Doub beta(const Doub z, const Doub w);

// End gamma.h.
// Begin contents of deviates.h
struct Expondev : Ran {
	Doub beta;
	Expondev(Doub bbeta, Ullong i);
	Doub dev();
};

struct Logisticdev : Ran {
	Doub mu,sig;
	Logisticdev(Doub mmu, Doub ssig, Ullong i);
	Doub dev();
};

struct Normaldev_BM : Ran {
	Doub mu,sig;
	Doub storedval;
	Normaldev_BM(Doub mmu, Doub ssig, Ullong i);
	Doub dev();
};

struct Cauchydev : Ran {
	Doub mu,sig;
	Cauchydev(Doub mmu, Doub ssig, Ullong i);
	Doub dev();
};

struct Normaldev : Ran {
	Doub mu,sig;
	Normaldev(Doub mmu, Doub ssig, Ullong i);
	Doub dev();
};

struct Gammadev : Normaldev {
	Doub alph, oalph, bet;
	Doub a1,a2;
	Gammadev(Doub aalph, Doub bbet, Ullong i);
	Doub dev();
};

struct Poissondev : Ran {
	Doub lambda, sqlam, loglam, lamexp, lambold;
	VecDoub logfact;
	Poissondev(Doub llambda, Ullong i);
	Int dev();
	Int dev(Doub llambda);
};

struct Binomialdev : Ran {
	Doub pp,p,pb,expnp,np,glnp,plog,pclog,sq;
	Int n,swch;
	Ullong uz,uo,unfin,diff,rltp;
	Int pbits[5];
	Doub cdf[64];
	Doub logfact[1024];
	Binomialdev(Int nn, Doub ppp, Ullong i);	
	Int dev();
};
// End deviates.h

}; // END OF NR3 NAMESPACE

// all ncclib entities are part of the ncc namespace
namespace ncc {

// Functions

void whoami();

#define ncc__error(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); exit(1);}

#define ncc__warning(message) \
{printf("WARNING: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__);}

void GetStrPropertyFromINIFile(const char *group, const char *prop, const char *defval, char *buf, const unsigned int bufsize, const char *filename);

int GetIntPropertyFromINIFile(const char *group, const char *prop, const int defval, const char *filename);

#ifdef _MSC_VER // this function uses the Microsoft specific strtok_s
float** QDReadMatrixFromFile(FILE* fp, unsigned int headerlines=0);
#endif // _MSC_VER block

#ifdef _MSC_VER // not sure why this doesn't work with gcc
vector<vector<float>*> *RSReadMatrixFromFile(char* filename, unsigned int headerlines=0);
#endif // cl only function

bool load(const char* filename, nr3::MatDoub *data, unsigned int headerlines=0);

bool logEntry(const char* filename, const char* entry, bool newlog=false);

struct vec3f {
// A rudimentary public access Cartesian 3-vector class.
	float x,y,z; // Cartesian components

	// constructors
	vec3f();
	vec3f(const float v);
	vec3f(const float vx, const float vy, const float vz);

	// operators
	vec3f operator +  (vec3f rhs);
	vec3f operator += (vec3f rhs);
	vec3f operator -  (vec3f rhs);
	vec3f operator -= (vec3f rhs);
	vec3f operator *  (float rhs);
	vec3f operator *= (float rhs);
	vec3f operator /  (float rhs);
	vec3f operator /= (float rhs);
	friend vec3f operator * (float a, vec3f v);

}; // END VEC3F
 
struct vec3d {
// A rudimentary public access Cartesian 3-vector class (double).
	double x,y,z; // Cartesian components

	// constructors
	vec3d();
	vec3d(const float v);
	vec3d(const float vx, const float vy, const float vz);

	// operators
	vec3d operator +  (vec3d rhs);
	vec3d operator += (vec3d rhs);
	vec3d operator -  (vec3d rhs);
	vec3d operator -= (vec3d rhs);
	vec3d operator *  (float rhs);
	vec3d operator *= (float rhs);
	vec3d operator /  (float rhs);
	vec3d operator /= (float rhs);
	friend vec3d operator * (float a, vec3d v);

}; // END VEC3D

// A namespace for RGB triplets with names, based on the X11 color specification
namespace rgb
{
	typedef unsigned char rgbTrip [3];
	// Red colors
	static const rgbTrip rIndianRed		=	{205, 92, 92};
	static const rgbTrip rLightCoral	=	{240,128,128};
	static const rgbTrip rSalmon		=	{250,128,114};
	static const rgbTrip rDarkSalmon	=	{233,150,122};
	static const rgbTrip rLightSalmon	=	{255,160,122};
	static const rgbTrip rRed		=	{255,  0,  0};
	static const rgbTrip rCrimson		=	{220, 20, 60};
	static const rgbTrip rFireBrick		=	{178, 34, 34};
	static const rgbTrip rDarkRed		=	{139,  0,  0};

	// Pink colors
	static const rgbTrip piPink		=	{255,192,203};
	static const rgbTrip piLightPink	=	{255,182,193};
	static const rgbTrip piHotPink		=	{255,105,180};
	static const rgbTrip piDeepPink		=	{255, 20,147};
	static const rgbTrip piMediumVioletRed	=	{199, 21,133};
	static const rgbTrip piPaleVioletRed	=	{219,112,147};

	// Orange colors
	static const rgbTrip oLightSalmon	=	{255,160,122};
	static const rgbTrip oCoral		=	{255,127, 80};
	static const rgbTrip oTomato		=	{255, 99, 71};
	static const rgbTrip oOrangeRed		=	{255, 69,  0};
	static const rgbTrip oDarkOrange	=	{255,140,  0};
	static const rgbTrip oOrange		=	{255,165,  0};

	// Yellow colors
	static const rgbTrip yGold		=	{255,215,  0};
	static const rgbTrip yYellow		=	{255,255,  0};
	static const rgbTrip yLightYellow	=	{255,255,224};
	static const rgbTrip yLemonChiffon	=	{255,250,205};
	static const rgbTrip yLightGoldenrod	=	{255,250,210};
	static const rgbTrip yPapayaWhip	=	{255,239,213};
	static const rgbTrip yMoccasin		=	{255,228,181};
	static const rgbTrip yPeachPuff		=	{255,218,185};
	static const rgbTrip yPaleGoldenrod	=	{238,232,170};
	static const rgbTrip yKhaki		=	{240,230,140};
	static const rgbTrip yDarkKhaki		=	{189,183,107};

	// Purple colors
	static const rgbTrip puLavender		=	{230,230,250};
	static const rgbTrip puThistle		=	{216,192,216};
	static const rgbTrip puPlum		=	{221,160,221};
	static const rgbTrip puViolet		=	{238,130,238};
	static const rgbTrip puOrchid		=	{218,112,214};
	static const rgbTrip puFuchsia		=	{255,  0,255};
	static const rgbTrip puMagenta		=	{255,  0,255};
	static const rgbTrip puMediumOrchid	=	{186, 85,211};
	static const rgbTrip puMediumPurple	=	{147,112,219};
	static const rgbTrip puBlueViolet	=	{138, 43,226};
	static const rgbTrip puDarkViolet	=	{148,  0,211};
	static const rgbTrip puDarkOrchid	=	{153, 50,204};
	static const rgbTrip puDarkMagenta	=	{139,  0,139};
	static const rgbTrip puPurple		=	{128,  0,128};
	static const rgbTrip puIndigo		=	{ 75,  0,130};
	static const rgbTrip puDarkSlate	=	{ 72, 61,139};
	static const rgbTrip puSlate		=	{106, 90,205};
	static const rgbTrip puMediumSlate	=	{123,104,238};

	// Green colors
	static const rgbTrip gGreenYellow	=	{173,255, 47};
	static const rgbTrip gChartreuse	=	{127,255,  0};
	static const rgbTrip gLawnGreen		=	{124,252,  0};
	static const rgbTrip gLime		=	{  0,255,  0};
	static const rgbTrip gLimeGreen		=	{ 50,255, 50};
	static const rgbTrip gPaleGreen		=	{152,251,152};
	static const rgbTrip gLightGreen	=	{144,238,144};
	static const rgbTrip gMediumSpring	=	{  0,250,154};
	static const rgbTrip gSpring		=	{  0,255,127};
	static const rgbTrip gMediumSea		=	{ 60,179,113};
	static const rgbTrip gSeaGreen		=	{ 46,139, 87};
	static const rgbTrip gForestGreen	=	{ 34,139, 34};
	static const rgbTrip gGreen		=	{  0,128,  0};
	static const rgbTrip gDarkGreen		=	{  0,100,  0};
	static const rgbTrip gYellowGreen	=	{154,205, 50};
	static const rgbTrip gOliveDrab		=	{107,142, 35};
	static const rgbTrip gOlive		=	{128,128,  0};
	static const rgbTrip gDarkOlive		=	{ 85,107, 47};
	static const rgbTrip gMediumAquamarine	=	{102,205,170};
	static const rgbTrip gDarkSeaGreen	=	{143,188,143};
	static const rgbTrip gLightSeaGreen	=	{ 32,178,170};
	static const rgbTrip gDarkCyan		=	{  0,139,139};
	static const rgbTrip gTeal		=	{  0,128,128};

	// Blue colors
	static const rgbTrip bAqua		=	{  0,255,255};
	static const rgbTrip bCyan		=	{  0,255,255};
	static const rgbTrip bLightCyan		=	{224,255,255};
	static const rgbTrip bPaleTurquoise	=	{175,238,238};
	static const rgbTrip bAquaMarine	=	{127,255,212};
	static const rgbTrip bTurquoise		=	{ 64,224,208};
	static const rgbTrip bMediumTurquoise	=	{ 72,209,204};
	static const rgbTrip bDarkTurqoise	=	{  0,206,209};
	static const rgbTrip bCadetBlue		=	{ 95,158,160};
	static const rgbTrip bSteelBlue		=	{ 70,130,180};
	static const rgbTrip bLightSteel	=	{176,196,222};
	static const rgbTrip bPowderBlue	=	{176,224,230};
	static const rgbTrip bLightBlue		=	{173,216,230};
	static const rgbTrip bSkyBlue		=	{135,206,235};
	static const rgbTrip bLightSkyBlue	=	{135,206,250};
	static const rgbTrip bDeepSkyBlue	=	{  0,191,255};
	static const rgbTrip bDodgerBlue	=	{ 30,144,255};
	static const rgbTrip bCornflourBlue	=	{100,149,237};
	static const rgbTrip bRoyalBlue		=	{ 65,105,225};
	static const rgbTrip bBlue		=	{  0,  0,255};
	static const rgbTrip bMediumBlue	=	{  0,  0,205};
	static const rgbTrip bNavy		=	{  0,  0,128};
	static const rgbTrip bMidnightBlue	=	{ 25, 25,112};

	// Brown colors
	static const rgbTrip brCornsilk		=	{255,248,220};
	static const rgbTrip brBlanchedAlmond	=	{255,235,205};
	static const rgbTrip brBisque		=	{255,228,196};
	static const rgbTrip brNavajoWhite	=	{255,222,173};
	static const rgbTrip brWheat		=	{245,222,179};
	static const rgbTrip brBurlyWood	=	{222,184,135};
	static const rgbTrip brTan		=	{210,180,140};
	static const rgbTrip brRosyBrown	=	{188,143,143};
	static const rgbTrip brSandyBrown	=	{244,164, 96};
	static const rgbTrip brGoldenrod	=	{218,165, 32};
	static const rgbTrip brDarkGoldenrod	=	{184,134, 11};
	static const rgbTrip brPeru		=	{205,133, 63};
	static const rgbTrip brChocolate	=	{210,105, 30};
	static const rgbTrip brSaddleBrown	=	{139, 69, 19};
	static const rgbTrip brSienna		=	{160, 82, 45};
	static const rgbTrip brBrown		=	{165, 42, 42};
	static const rgbTrip brMaroon		=	{128,  0,  0};

	// White colors
	static const rgbTrip wWhite		=	{255,255,255};
	static const rgbTrip wSnow		=	{255,250,250};
	static const rgbTrip wHoneydew		=	{240,255,240};
	static const rgbTrip wMintCream		=	{245,255,250};
	static const rgbTrip wAzure		=	{240,255,255};
	static const rgbTrip wAliceBlue		=	{240,248,255};
	static const rgbTrip wGhostWhite	=	{248,248,255};
	static const rgbTrip wWhiteSmoke	=	{245,245,245};
	static const rgbTrip wSeashell		=	{255,245,238};
	static const rgbTrip wBeige		=	{245,245,220};
	static const rgbTrip wOldLace		=	{253,245,230};
	static const rgbTrip wFloralWhite	=	{255,250,240};
	static const rgbTrip wIvory		=	{255,255,240};
	static const rgbTrip wAntiqueWhite	=	{250,235,215};
	static const rgbTrip wLinen		=	{250,240,230};
	static const rgbTrip wLavenderBlush	=	{255,240,245};
	static const rgbTrip wMistyRose		=	{255,228,225};

	// Gray colors
	static const rgbTrip grGainsboro	=	{220,220,220};
	static const rgbTrip grLightGray	=	{211,211,211};
	static const rgbTrip grSilver		=	{192,192,192};
	static const rgbTrip grDarkGray		=	{169,169,169};
	static const rgbTrip grGray		=	{128,128,128};
	static const rgbTrip grDimGray		=	{105,105,105};
	static const rgbTrip grLightSlateGray	=	{119,136,153};
	static const rgbTrip grSlateGray	=	{112,128,144};
	static const rgbTrip grDarkSlateGray	=	{ 47, 79, 79};
	static const rgbTrip grBlack		=	{  0,  0,  0};
}; // END RGB NAMESPACE

}; // END THE NCC NAMESPACE

#endif // ncclib include guard
