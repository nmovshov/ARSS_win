////////////////////////////////////////////////////////////////////////////////
//NCCLIB	Naor's C collection library
//
//		This is the implementation.
// Contents:
// 	0) #include ncclib.h. Forward declarations, platform specific macros, and
// 	   system #includes, e.g. windows.h and Microsoft compiler macros from the
// 	   CRT.
//
// 	1) Modified nr3.h. This is a modified version of the Numerical Recipes
// 	   third edition main dependency. To simplify transportability, I decided
// 	   to copy nr3.h and other useful nr routines into ncclib.cpp, rather than
// 	   #include their respective .h files. Since nr doesn't bother with forward
// 	   declarations I had to create those for all nr3 routines and put them
// 	   in ncclib.h, then make the original contents of nr3.h into definitions.
// 	   The contents of nr3.h are surrounded in a nr3 namespace.
//
// 	  1a) The original nr3.h, inside an nr3 namespace. These include a few
// 	      useful macro names for inline functions (MAX,MIN, etc.), an error
// 	      macro (throw), a quiet NaN variable, and an efficient implementation
// 	      of templated vector and matrix classes.
// 	  1b) The contents of ran.h, also inside the nr3 namespace. Much better than
// 	      the built-in rand(). The usage is:
// 		nr3::Ran ran(seed);
// 		ran.int64(); for a random uint64
// 		ran.doub();  for a random double in (0,1)
// 		1 + ran.int64() % (n-1) for a random uint64 between 1 and n
// 	  1c) The contents of gamma.h
// 	  1d) The contents of deviates.h. Some useful probability distributions.
// 	      Usage:
// 		nr3::Normaldev norm(mu.sig,seed);
// 		norm.dev(); for a normal deviate (copied from deviates.h)
//
// 	2) WHOMAI(). An expanded "hello world" message with environment info.
//
// 	3) NCC__ERROR() and NCC__WARNING(). Simple macro printf and exit(1); header only!
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
// 	 5) VEC3F and VEC3D. Rudimentary public access Cartesian 3-vector class.
//
// 	 6) RGB. A namespace for named colors based on the X11 color spec. Header only!
// 	
// Author: Me (Naor Movshovitz)
////////////////////////////////////////////////////////////////////////////////

#include "ncclib.h"

// To request lint level warnings on Microsoft compilers, uncomment the following
// #ifdef and the corresponding block at the end of the file, and #define the LINT macro.
/*#ifdef LINT
#pragma warning(push,4)
#pragma warning(disable:4100) // unreferenced formal parameter
#endif*/

// All NR3 entities are part of the nr3 namespace
namespace nr3 {
// Floating Point Exceptions for Microsoft compilers
#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp( 0, 0 );
		cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
		_controlfp( cw, MCW_EM );
	}
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */


// End of original nr3.h
// Begin contents of ran.h
Doub Ran::doub() { return 5.42101086242752217E-20 * int64(); }

Uint Ran::int32() { return (Uint)int64(); }

Ullong Ran::int64() {
	u = u * 2862933555777941757LL + 7046029254386353087LL;
	v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
	w = 4294957665U*(w & 0xffffffff) + (w >> 32);
	Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
	return (x + v) ^ w;
}

Ran::Ran( Ullong j ) : v(4101842887655102017LL), w(1)
{
	u = j ^ v; int64();
	v = u; int64();
	w = v; int64();
}

Ranq1::Ranq1( Ullong j ) : v(4101842887655102017LL)
{
	v ^= j;
	v = int64();
}

Ullong Ranq1::int64()
{
	v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
	return v * 2685821657736338717LL;
}

Doub Ranq1::doub() { return 5.42101086242752217E-20 * int64(); }

Uint Ranq1::int32() { return (Uint)int64(); }

Ranq2::Ranq2( Ullong j ) : v(4101842887655102017LL), w(1)
{
	v ^= j;
	w = int64();
	v = int64();
}

Ullong Ranq2::int64()
{
	v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
	w = 4294957665U*(w & 0xffffffff) + (w >> 32);
	return v ^ w;
}

Doub Ranq2::doub() { return 5.42101086242752217E-20 * int64(); }

Uint Ranq2::int32() { return (Uint)int64(); }

Ullong Ranhash::int64( Ullong u )
{
	Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
	v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
	v *= 4768777513237032717LL;
	v ^= v << 20; v ^= v >> 41; v ^= v << 5;
	return  v;
}

Uint Ranhash::int32( Ullong u ) { return (Uint)(int64(u) & 0xffffffff) ; }

Doub Ranhash::doub( Ullong u ) { return 5.42101086242752217E-20 * int64(u); }

Ranbyte::Ranbyte( Int u )
{
	v = 2244614371U ^ u;
	for (i=0; i<256; i++) {s[i] = i;}
	for (j=0, i=0; i<256; i++) {
		ss = s[i];
		j = (j + ss + (v >> 24)) & 0xff;
		s[i] = s[j]; s[j] = ss;
		v = (v << 24) | (v >> 8);
	}
	i = j = 0;
	for (Int k=0; k<256; k++) int8();
}

unsigned char Ranbyte::int8()
{
	i = (i+1) & 0xff;
	ss = s[i];
	j = (j+ss) & 0xff;
	s[i] = s[j]; s[j] = ss;
	return (unsigned char)(s[(s[i]+s[j]) & 0xff]);
}

Uint Ranbyte::int32()
{
	v = 0;
	for (int k=0; k<4; k++) {
		i = (i+1) & 0xff;
		ss = s[i];
		j = (j+ss) & 0xff;
		s[i] = s[j]; s[j] = ss;
		v = (v << 8) | s[(s[i]+s[j]) & 0xff];
	}
	return v;
}

Doub Ranbyte::doub()
{
	return 2.32830643653869629E-10 * ( int32() +
		2.32830643653869629E-10 * int32() );
}

Ranfib::Ranfib( Ullong j ) : inext(0), inextp(31)
{
	Ranq1 init(j);
	for (int k=0; k<55; k++) dtab[k] = init.doub();
}

Doub Ranfib::doub()
{
	if (++inext == 55) inext = 0;
	if (++inextp == 55) inextp = 0;
	dd = dtab[inext] - dtab[inextp];
	if (dd < 0) dd += 1.0;
	return (dtab[inext] = dd);
}

unsigned long Ranfib::int32() { return (unsigned long)(doub() * 4294967295.0); }

Ranlim32::Ranlim32( Uint j ) : v(2244614371U), w1(521288629U), w2(362436069U)
{
	u = j ^ v; int32();
	v = u; int32();
}

Uint Ranlim32::int32()
{
	u = u * 2891336453U + 1640531513U;
	v ^= v >> 13; v ^= v << 17; v ^= v >> 5;
	w1 = 33378 * (w1 & 0xffff) + (w1 >> 16);
	w2 = 57225 * (w2 & 0xffff) + (w2 >> 16);
	Uint x = u ^ (u << 9); x ^= x >> 17; x ^= x << 6;
	Uint y = w1 ^ (w1 << 17); y ^= y >> 15; y ^= y << 5;
	return (x + v) ^ (y + w2);
}

Doub Ranlim32::doub() { return 2.32830643653869629E-10 * int32(); }

Doub Ranlim32::truedoub()
{
	return 2.32830643653869629E-10 * ( int32() +
		2.32830643653869629E-10 * int32() );
}

// End of Ran.h.
// Begin contents of gamma.h
Doub gammln(const Doub xx) {
	Int j;
	Doub x,tmp,y,ser;
	static const Doub cof[14]={57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}

Doub factrl(const Int n) {
	static VecDoub a(171);
	static Bool init=true;
	if (init) {
		init = false;
		a[0] = 1.;
		for (Int i=1;i<171;i++) a[i] = i*a[i-1];
	}
	if (n < 0 || n > 170) throw("factrl out of range");
	return a[n];
}

Doub factln(const Int n) {
	static const Int NTOP=2000;
	static VecDoub a(NTOP);
	static Bool init=true;
	if (init) {
		init = false;
		for (Int i=0;i<NTOP;i++) a[i] = gammln(i+1.);
	}
	if (n < 0) throw("negative arg in factln");
	if (n < NTOP) return a[n];
	return gammln(n+1.);
}

Doub bico(const Int n, const Int k) {
	if (n<0 || k<0 || k>n) throw("bad args in bico");
	if (n<171) return floor(0.5+factrl(n)/(factrl(k)*factrl(n-k)));
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

Doub beta(const Doub z, const Doub w) {
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}

// End gamma.h.
// Begin contents of deviates.h
Expondev::Expondev( Doub bbeta, Ullong i ) : Ran(i), beta(bbeta) { }

Doub Expondev::dev()
{
	Doub u;
	do u = doub(); while (u == 0.);
	return -log(u)/beta;
}

Logisticdev::Logisticdev( Doub mmu, Doub ssig, Ullong i ) : Ran(i), mu(mmu), sig(ssig) { }

Doub Logisticdev::dev()
{
	Doub u;
	do u = doub(); while (u*(1.-u) == 0.);
	return mu + 0.551328895421792050*sig*log(u/(1.-u));
}

Normaldev_BM::Normaldev_BM( Doub mmu, Doub ssig, Ullong i ) : Ran(i), mu(mmu), sig(ssig), storedval(0.) { }

Doub Normaldev_BM::dev()
{
	Doub v1,v2,rsq,fac;
	if (storedval == 0.) {
		do {
			v1=2.0*doub()-1.0;
			v2=2.0*doub()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		storedval = v1*fac;
		return mu + sig*v2*fac;
	} else {
		fac = storedval;
		storedval = 0.;
		return mu + sig*fac;
	}
}

Cauchydev::Cauchydev( Doub mmu, Doub ssig, Ullong i ) : Ran(i), mu(mmu), sig(ssig) { }

Doub Cauchydev::dev()
{
	Doub v1,v2;
	do {
		v1=2.0*doub()-1.0;
		v2=doub();
	} while (SQR(v1)+SQR(v2) >= 1. || v2 == 0.);
	return mu + sig*v1/v2;
}

Normaldev::Normaldev( Doub mmu, Doub ssig, Ullong i ) : Ran(i), mu(mmu), sig(ssig) { }

Doub Normaldev::dev()
{
	Doub u,v,x,y,q;
	do {
		u = doub();
		v = 1.7156*(doub()-0.5);
		x = u - 0.449871;
		y = abs(v) + 0.386595;
		q = SQR(x) + y*(0.19600*y-0.25472*x);
	} while (q > 0.27597
		&& (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
	return mu + sig*v/u;
}

Gammadev::Gammadev( Doub aalph, Doub bbet, Ullong i ) : Normaldev(0.,1.,i), alph(aalph), oalph(aalph), bet(bbet)
{
	if (alph <= 0.) throw("bad alph in Gammadev");
	if (alph < 1.) alph += 1.;
	a1 = alph-1./3.;
	a2 = 1./sqrt(9.*a1);
}

Doub Gammadev::dev()
{
	Doub u,v,x;
	do {
		do {
			x = Normaldev::dev();
			v = 1. + a2*x;
		} while (v <= 0.);
		v = v*v*v;
		u = doub();
	} while (u > 1. - 0.331*SQR(SQR(x)) &&
		log(u) > 0.5*SQR(x) + a1*(1.-v+log(v)));
	if (alph == oalph) return a1*v/bet;
	else {
		do u=doub(); while (u == 0.);
		return pow(u,1./oalph)*a1*v/bet;
	}
}

Poissondev::Poissondev( Doub llambda, Ullong i ) : Ran(i), lambda(llambda), logfact(1024,-1.), lambold(-1.) { }

Int Poissondev::dev()
{
	Doub u,u2,v,v2,p,t,lfac;
	Int k;
	if (lambda < 5.) {
		if (lambda != lambold) lamexp=exp(-lambda);
		k = -1;
		t=1.;
		do {
			++k;
			t *= doub();
		} while (t > lamexp);
	} else {
		if (lambda != lambold) {
			sqlam = sqrt(lambda);
			loglam = log(lambda);
		}
		for (;;) {
			u = 0.64*doub();
			v = -0.68 + 1.28*doub();
			if (lambda > 13.5) {
				v2 = SQR(v);
				if (v >= 0.) {if (v2 > 6.5*u*(0.64-u)*(u+0.2)) continue;}
				else {if (v2 > 9.6*u*(0.66-u)*(u+0.07)) continue;}
			}
			k = Int(floor(sqlam*(v/u)+lambda+0.5));
			if (k < 0) continue;
			u2 = SQR(u);
			if (lambda > 13.5) {
				if (v >= 0.) {if (v2 < 15.2*u2*(0.61-u)*(0.8-u)) break;}
				else {if (v2 < 6.76*u2*(0.62-u)*(1.4-u)) break;}
			}
			if (k < 1024) {
				if (logfact[k] < 0.) logfact[k] = gammln(k+1.);
				lfac = logfact[k];
			} else lfac = gammln(k+1.);
			p = sqlam*exp(-lambda + k*loglam - lfac);
			if (u2 < p) break;
		}
	}
	lambold = lambda;
	return k;
}

Int Poissondev::dev( Doub llambda )
{
	lambda = llambda;
	return dev();
}

Binomialdev::Binomialdev( Int nn, Doub ppp, Ullong i ) : Ran(i), pp(ppp), n(nn)
{
	Int j;
	pb = p = (pp <= 0.5 ? pp : 1.0-pp);
	if (n <= 64) {
		uz=0;
		uo=0xffffffffffffffffLL;
		rltp = 0;
		for (j=0;j<5;j++) pbits[j] = 1 & ((Int)(pb *= 2.));
		pb -= floor(pb);
		swch = 0;
	} else if (n*p < 30.) {
		cdf[0] = exp(n*log(1-p));
		for (j=1;j<64;j++) cdf[j] =  cdf[j-1] + exp(gammln(n+1.)
			-gammln(j+1.)-gammln(n-j+1.)+j*log(p)+(n-j)*log(1.-p));
			swch = 1;
	} else {
		np = n*p;
		glnp=gammln(n+1.);
		plog=log(p);
		pclog=log(1.-p);
		sq=sqrt(np*(1.-p));
		if (n < 1024) for (j=0;j<=n;j++) logfact[j] = gammln(j+1.);
		swch = 2;
	}
}

Int Binomialdev::dev()
{
	Int j,k,kl,km;
	Doub y,u,v,u2,v2,b;
	if (swch == 0) {
		unfin = uo;
		for (j=0;j<5;j++) {
			diff = unfin & (int64()^(pbits[j]? uo : uz));
			if (pbits[j]) rltp |= diff;
			else rltp = rltp & ~diff;
			unfin = unfin & ~diff;
		}
		k=0;
		for (j=0;j<n;j++) {
			if (unfin & 1) {if (doub() < pb) ++k;}
			else {if (rltp & 1) ++k;}
			unfin >>= 1;
			rltp >>= 1;
		}
	} else if (swch == 1) {
		y = doub();
		kl = -1;
		k = 64;
		while (k-kl>1) {
			km = (kl+k)/2;
			if (y < cdf[km]) k = km;
			else kl = km;
		}
	} else {
		for (;;) {
			u = 0.645*doub();
			v = -0.63 + 1.25*doub();
			v2 = SQR(v);
			if (v >= 0.) {if (v2 > 6.5*u*(0.645-u)*(u+0.2)) continue;}
			else {if (v2 > 8.4*u*(0.645-u)*(u+0.1)) continue;}
			k = Int(floor(sq*(v/u)+np+0.5));
			if (k < 0) continue;
			u2 = SQR(u);
			if (v >= 0.) {if (v2 < 12.25*u2*(0.615-u)*(0.92-u)) break;}
			else {if (v2 < 7.84*u2*(0.615-u)*(1.2-u)) break;}
			b = sq*exp(glnp+k*plog+(n-k)*pclog
				- (n < 1024 ? logfact[k]+logfact[n-k]
			: gammln(k+1.)+gammln(n-k+1.)));
			if (u2 < b) break;
		}
	}
	if (p != pp) k = n - k;
	return k;
}
// End deviates.h
}; // End the NR3 NAMESPACE

// all ncclib entities are part of the ncc namespace
namespace ncc {

// Functions

void whoami()
{
//WHOAMI  Get some info on the current machine
//  WHOAMI prints out a little "hello world" message with some additional
//  information on the running environment.

	const char *host=NULL;
	const char *os=NULL;
	const char *user=NULL;
	const char *endiandness=NULL;
	const char *macharch=NULL;

#if defined _WIN32
	host     = getenv("computername");
	os       = getenv("OS");
	user     = getenv("username");
	macharch = getenv("processor_architecture");
#elif defined __linux__ or defined __unix__
	host     = getenv("HOST");
	os       = getenv("OSTYPE");
	user     = getenv("USER");
	macharch = getenv("MACHTYPE");
#elif defined __APPLE__
	host     = getenv("HOST");
	os       = getenv("OSTYPE");
	user     = getenv("USER");
	macharch = getenv("MACHTYPE");
#else
	host     = "unknown";
	user     = "the user";
	os       = "unknown";
	macharch = "unknown";
#endif

#ifdef __GNUC__
#if   __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
	endiandness="little-endian";
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
	endiandness="big-endian";
#elif __BYTE_ORDER__ == __ORDER_PDP_ENDIAN__
	endiandness="mixed-endian";
#else
	endiandness="unknown-endian";
#endif
#else
	endiandness="unknown-endian";
#endif

	printf("Alo world, and a special alo to %s.\n",user);
	printf("My name is %s, and I am a %s machine running %s.\n",host,macharch,os);
	printf("I am currently compiling with %s %u-bits.\n",endiandness,8*sizeof(void*));

} // END FUNCTION WHOAMI

void GetStrPropertyFromINIFile(const char *group, const char *prop, const char *defval, char *buf, const unsigned int bufsize, const char *filename)
{
// This function mimics the behavior of GetPrivateProfileStr() from the Windows
// API. The function looks inside the .ini file (old-school windows input
// parameters) for the specified property name, inside the specified group name.
// This provides some flexibility in defining groups of related parameters with
// short names. The Windows function works well. Since this function is typically
// only called at the beginning of the program, efficiency is not critical and the
// file will be opened and closed every time this function is called. Reasonable
// input checking is included, but a dick user can easily break through, so don't
// be a dick. Also, on non-windows platforms this function is case sensitive.
// Oh and in case you forgot, the format of a .ini file is:
// ; comment
// [group]
// property=value
// ...
// and the winapi function takes a full path to work!

#ifdef _WIN32 // might as well use the original is it's there...
	GetPrivateProfileString(group, prop, defval, buf, bufsize, filename);
#else
	char line[256];
	char *tok;
	strcpy(buf,defval);

	// first, verify access to the file
	ifstream fid(filename);
	if (!fid.good()) return; // note: no error or warning

	// look for group
	fid.seekg(0,ios::beg);
	char gs[256];
	sprintf(gs,"%s%s%s","[",group,"]");
	while (fid.good())
	{
		fid.getline(line,256);
		if (strcmp(line,gs)==0) break;
	}
	if (fid.eof()) // group not found (or is empty)
	{
		fid.close();
		return;
	}

	// look inside the group
	while (fid.good())
	{
		fid.getline(line,256);
		if (strchr(line,'[')) // fell into next group...
		{
			fid.close();
			return;
		}
		tok = strtok(line," \t =");
		if (tok == NULL) continue;
		if (strcmp(tok,prop)!=0) continue;
		char *eqpos=strchr(line,'=');
		tok = strtok(eqpos,"= \t");
		if (tok == NULL) continue;
		strcpy(buf,tok);
		break;
	}

	fid.close();

#endif // windows/linux branching
} // END FUNCTION GETSTRPROPERTYFROMFILE

int GetIntPropertyFromINIFile(const char *group, const char *prop, const int defval, const char *filename)
{
// This function mimics the behavior of GetPrivateProfileInt() from the Windows
// API. The function looks inside the .ini file (old-school windows input
// parameters) for the specified property name, inside the specified group name.
// This provides some flexibility in defining groups of related parameters with
// short names. The Windows function works well. Since this function is typically
// only called at the beginning of the program, efficiency is not critical and the
// file will be opened and closed every time this function is called. Reasonable
// input checking is included, but a dick user can easily break through, so don't
// a dick. Also, on non-windows platforms this function is case sensitive.
// Oh and in case you forgot, the format of a .ini file is:
// ; comment
// [group]
// property=value
// ...
// and the winapi function takes a full path to work!

#ifdef _WIN32 // might as well use the original is it's there...
	return GetPrivateProfileInt(group, prop, defval, filename);
#else
	char line[256];
	char *tok;
	int intprop=defval;

	// first, verify access to the file
	ifstream fid(filename);
	if (!fid.good()) return intprop; // note: no error or warning

	// look for group
	fid.seekg(0,ios::beg);
	char gs[256];
	sprintf(gs,"%s%s%s","[",group,"]");
	while (fid.good())
	{
		fid.getline(line,256);
		if (strcmp(line,gs)==0) break;
	}
	if (fid.eof()) // group not found (or is empty)
	{
		fid.close();
		return intprop;
	}

	// look inside the group
	while (fid.good())
	{
		fid.getline(line,256);
		if (strchr(line,'[')) // fell into next group...
		{
			fid.close();
			return intprop;
		}
		tok = strtok(line," \t =");
		if (tok == NULL) continue;
		if (strcmp(tok,prop)!=0) continue;
		char *eqpos=strchr(line,'=');
		tok = strtok(eqpos,"= \t");
		if (tok == NULL) continue;
		intprop=(int)atof(tok);
		break;
	}

	fid.close();
	return intprop;

#endif // windows/linux branching
} // END FUNCTION GETINTPROPERTYFROMFILE

#ifdef _MSC_VER // this function uses the Microsoft specific strtok_s
float** QDReadMatrixFromFile(FILE* fp, unsigned int headerlines)
{
/*
A Quick and Dirty matrix reader using c-style arrays and fscans
This function takes a file pointer and returns a pointer to a c-style "matrix"
(float**) with values read from the file. An optional integer number of header
lines can be specified and they will be skipped. Note that this function is
quite brittle. Specifically, the user must ensure that fp is a valid pointer
to an open text file, and that, except for the header, the file contains nothing
but a rectangular array of real numbers.
The numbers will be read as floats even if they are ints. The number of columns
in the resulting matrix is determined from the number of records in the first
non-header line, and the number of rows is determined by the number of non-empty
lines in the file. Unfortunately, once the function exits this information is
lost. The caller can operate on the return array, but can't know its size. This
limitation makes the usefulness of this function doubtful. Mostly its here to
demonstrate the techniques of reading and parsing with old style c IO.
*/

	// First, skip header lines
	rewind(fp);
	char line[255];
	int count=headerlines;
	while (count--) fgets(line,255,fp);

	// Next, determine number of matrix columns, based on first data line.
	fgets(line,255,fp);
	int ncols=0;
	char* token=NULL;
	char* context=NULL;
	char seps[]=" ,\t\n"; // list of delimiters
	token=strtok_s(line,seps,&context);
	while (token!=NULL)	{ncols++; token=strtok_s(NULL,seps,&context);}

	// Next, determine number of matrix rows.
	int nrows=1;
	*line=NULL;
	while (!feof(fp)) 
	{
		fgets(line,255,fp);
		if (strcmp(line,"\n")) // don't count empty lines
			nrows++;
	}

	// Now, allocate memory for a matrix of nrows-by-ncols.
	float** matrix=NULL;
	matrix = new float*[nrows];
	for (int j=0; j<nrows; j++)
		matrix[j] = new float[ncols];

	// Go back to first data row,
	rewind(fp);
	count=headerlines;
	while (count--) fgets(line,255,fp);

	// and read matrix line by line.
	// At this point the function will fail with anything but a rectangular numeric array.
	for (int j=0; j<nrows; j++)
		for (int k=0; k<ncols; k++)
			fscanf(fp,"%f",matrix[j]+k);

	return matrix;
} // END OF QDREADMATRIXFROMFILE
#endif // _MSC_VER block

#ifdef _MSC_VER // not sure why this doesn't work with gcc
vector<vector<float>*> *RSReadMatrixFromFile(char* filename, unsigned int headerlines)
{
/*
A reasonably safe matrix reader using c++ STL
This function open a stream to the file filename and attempts to read from it
assuming it contains a rectangular array of numbers. If successful a matrix
(c++ vector of vectors) is returned. Whatever was the numerical type in the file,
the returned matrix is float. Some safeguards are attempted, but it's the user's
responsibility to make sure the file contains, apart from an optional header, nothing
but numbers in a rectangular array. The number of columns is determined by the
first non-header line. Otherwise a warning is shown and the function returns NULL.
*/

	// First, check the the input file is accessible.
	ifstream fp(filename);
	if (!fp.good()) {cout<<"error opening file"<<endl; return NULL;}

	// OK, since the file is here, start by skipping the header.
	char line[255];
	int count=headerlines;
	while (count--) fp.getline(line,255);

	// Next, determine number of matrix columns, based on first data line.
	// (Try to skip accidentally empties.)
	do fp.getline(line,255); while (strcmp(line,"\0")==0);
	int ncols=0;
	char* token=NULL;
	char* context=NULL;
	char seps[]=" ,\t\n"; // list of delimiters
	token=strtok_s(line,seps,&context);
	float x;
	while (token!=NULL && sscanf_s(token,"%f",&x))
	{
		ncols++;
	       	token=strtok_s(NULL,seps,&context);
	}

	// Next, determine number of matrix rows. (Non-empty lines in file.)
	int nrows=1;
	*line=NULL; // (strtok mangles line beyond recognition)
	while (fp.good()) 
	{
		fp.getline(line,255);
		if (strcmp(line,"\0")) // don't count empty lines
			nrows++;
	}

	// Now that we know the expected matrix size, we can preallocate memory
	// and initialized to zero.
	vector<vector<float>*> *matrix=NULL;
	matrix = new vector<vector<float>*>(nrows,NULL);
	if (matrix==NULL) {cout<<"error allocating memory"<<endl; return NULL;}
	for (int k=0; k<nrows; k++){
		(matrix->at(k)) = (new  vector<float>(ncols,0));
		if ((matrix->at(k))==NULL) {cout<<"error allocating memory"<<endl; return NULL;}
	}

	// Memory allocated, time to fill with numbers. Start back at first data line,
	fp.clear(); // crucial! since eof was reached the stream is in fail mode and seek would have no effect!
	fp.seekg(0,ios::beg);
	count=headerlines;
	while (count--) fp.getline(line,255);

	// and read and fill line by line.
	for (int j=0; j<nrows; j++){
		for (int k=0; k<ncols; k++){
			fp>>skipws>>(matrix->at(j)->at(k));
			if (fp.fail()) 
			{cout<<"error reading file - may not be rectangular or contain non-numeric data"<<endl;
			return NULL;}
		}
		fp.getline(line,255);
	}

	// That's it. MATRIX should be filled with values read from the file.
	fp.close();
	return matrix;
} // END RSREADMATRIXFROMFILE
#endif // cl only function

bool load(const char* filename, nr3::MatDoub *data, unsigned int headerlines)
{
/*
A reasonably safe matrix reader using nr3 matrix types.
This function open a stream to the file filename and attempts to read from it
assuming it contains a rectangular array of numbers. If successful a matrix
is written to data, which must be a valid pointer to an existing nr3::MatDoub.
Whatever was the numerical type in the file, the output matrix is double. Some
safeguards are attempted, but it's the user's responsibility to make sure the 
file contains, apart from an optional header, nothing but numbers in a rectangular
array. The number of columns is determined by the first non-header line. The return
parameter indicates success.
*/

	// First, check that the input file is accessible.
	ifstream fp(filename);
	if (!fp.good()) return false;

	// And that the buffer is ready to receive output
	if (data==NULL) return false;

	// OK, since the file is here, start by skipping the header.
	char line[255];
	int count=headerlines;
	while (count--) fp.getline(line,255);

	// Next, determine number of matrix columns, based on first data line.
	// (Try to skip accidentally empties.)
	do fp.getline(line,255); while (strcmp(line,"\0")==0);
	int ncols=0;
	char* token=NULL;
	char* context=NULL;
	char seps[]=" ,\t\n"; // list of delimiters
#ifdef _MSC_VER // use secure version with Microsoft's compiler
	token=strtok_s(line,seps,&context);
	float x;
	while (token!=NULL && sscanf_s(token,"%f",&x))
	{
		ncols++;
	       	token=strtok_s(NULL,seps,&context);
	}
#else // use non-secure version with gcc etc.
	token=strtok(line,seps);
	float x;
	while (token!=NULL && sscanf(token,"%f",&x))
	{
		ncols++;
	       	token=strtok(NULL,seps);
	}
#endif

	// Next, determine number of matrix rows. (Non-empty lines in file.)
	int nrows=1;
//	line="\0"; // (strtok mangles line beyond recognition)
	while (fp.good()) 
	{
		fp.getline(line,255);
		if (strcmp(line,"\0")) // don't count empty lines
			nrows++;
	}

	// Now that we know the expected matrix size, we can preallocate memory
	data->resize(nrows,ncols);

	// Memory allocated, time to fill with numbers. Start back at first data line,
	fp.clear(); // crucial! Since eof was reached the stream is in fail mode and seek would have no effect!
	fp.seekg(0,ios::beg);
	count=headerlines;
	while (count--) fp.getline(line,255);

	// and read and fill line by line.
	for (int j=0; j<nrows; j++){
		for (int k=0; k<ncols; k++){
			fp>>skipws>>((*data)[j][k]);
			if (fp.fail()) 
			{cout<<"error reading file - may not be rectangular or contain non-numeric data"<<endl;
			return false;}
		}
		fp.getline(line,255);
	}

	// That's it. *DATA should be filled with values read from the file.
	fp.close();
	
	return true;
} // END LOAD

bool logEntry(const char* filename, const char* entry, bool newlog)
{
	if (newlog) // make sure file doesn't already exist
	{
		ifstream logfile(filename);
		if (logfile.is_open()) return false; // file will be closed by ~ifstream
	}
	ofstream logfile(filename, ios::app);
	if (!logfile.is_open()) return false;
	logfile << entry << endl;
	logfile.close();

	return true;
} // LOGENTRY
// Classes

// VEC3F - A rudimentary public access Cartesian 3-vector class.
	// constructors
	vec3f::vec3f() {x=0.0f; y=0.0f; z=0.0f;}
	vec3f::vec3f(const float v) {x=y=z=v;}
	vec3f::vec3f(const float vx, const float vy, const float vz) {x=vx; y=vy; z=vz;}

	// operators
	vec3f vec3f::operator + (vec3f rhs) {
		vec3f tmp;
	       	tmp.x=x+rhs.x;
		tmp.y=y+rhs.y;
		tmp.z=z+rhs.z;
		return tmp;
	}

	vec3f vec3f::operator += (vec3f rhs) {
		x+=rhs.x;
		y+=rhs.y;
		z+=rhs.z;
		return *this;
	}

	vec3f vec3f::operator - (vec3f rhs) {
		vec3f tmp;
	       	tmp.x=x-rhs.x;
		tmp.y=y-rhs.y;
		tmp.z=z-rhs.z;
		return tmp;
	}

	vec3f vec3f::operator -= (vec3f rhs) {
		x-=rhs.x;
		y-=rhs.y;
		z-=rhs.z;
		return *this;
	}

	vec3f vec3f::operator * (float rhs) {
		vec3f tmp;
		tmp.x=x*rhs;
		tmp.y=y*rhs;
		tmp.z=z*rhs;
		return tmp;
	}

	vec3f vec3f::operator *= (float rhs) {
		x*=rhs;
		y*=rhs;
		z*=rhs;
		return *this;
	}

	vec3f vec3f::operator / (float rhs) {
		vec3f tmp;
		tmp.x=x/rhs;
		tmp.y=y/rhs;
		tmp.z=z/rhs;
		return tmp;
	}

	vec3f vec3f::operator /= (float rhs) {
		x/=rhs;
		y/=rhs;
		z/=rhs;
		return *this;
	}

	vec3f operator * (float a, vec3f v) {
		v.x*=a;
		v.y*=a;
		v.z*=a;
		return v;
	}

 // END VEC3F
 
// VEC3D - A rudimentary public access Cartesian 3-vector class (double).
	// constructors
	vec3d::vec3d() {x=0.0f; y=0.0f; z=0.0f;}
	vec3d::vec3d(const float v) {x=y=z=v;}
	vec3d::vec3d(const float vx, const float vy, const float vz) {x=vx; y=vy; z=vz;}

	// operators
	vec3d vec3d::operator + (vec3d rhs) {
		vec3d tmp;
	       	tmp.x=x+rhs.x;
		tmp.y=y+rhs.y;
		tmp.z=z+rhs.z;
		return tmp;
	}

	vec3d vec3d::operator += (vec3d rhs) {
		x+=rhs.x;
		y+=rhs.y;
		z+=rhs.z;
		return *this;
	}

	vec3d vec3d::operator - (vec3d rhs) {
		vec3d tmp;
	       	tmp.x=x-rhs.x;
		tmp.y=y-rhs.y;
		tmp.z=z-rhs.z;
		return tmp;
	}

	vec3d vec3d::operator -= (vec3d rhs) {
		x-=rhs.x;
		y-=rhs.y;
		z-=rhs.z;
		return *this;
	}

	vec3d vec3d::operator * (float rhs) {
		vec3d tmp;
		tmp.x=x*rhs;
		tmp.y=y*rhs;
		tmp.z=z*rhs;
		return tmp;
	}

	vec3d vec3d::operator *= (float rhs) {
		x*=rhs;
		y*=rhs;
		z*=rhs;
		return *this;
	}

	vec3d vec3d::operator / (float rhs) {
		vec3d tmp;
		tmp.x=x/rhs;
		tmp.y=y/rhs;
		tmp.z=z/rhs;
		return tmp;
	}

	vec3d vec3d::operator /= (float rhs) {
		x/=rhs;
		y/=rhs;
		z/=rhs;
		return *this;
	}

	vec3d operator * (float a, vec3d v) {
		v.x*=a;
		v.y*=a;
		v.z*=a;
		return v;
	}

// END VEC3D

}; // END THE NCC NAMESPACE

// End lint level warnings
/*#ifdef LINT
#pragma warning(pop)
#endif*/
