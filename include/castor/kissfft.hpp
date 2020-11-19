/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : kissfft.hpp                                   |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 31.10.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Original kissfft library in one C++ .hpp file |
 |  `---'  |                for fourrier transform of castor full matrix. |
 +========================================================================+
 */

/*
 *  This fully integrated files are part and property of KISSFFT library :
 *   - "kiss_fft_log.h"
 *   - "kiss_fft.h"
 *   - "kiss_fftr.h"
 *   - "_kiss_fft_guts.h"
 *   - "kiss_fftr.c"
 *   - "kiss_fft.c"
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  SPDX-License-Identifier: BSD-3-Clause.
 *  For more informations :
 *  https://github.com/mborgerding/kissfft
 */

#pragma once
#define CASTOR_KISSFFT_HPP
#include <castor/matrix.hpp>


//==========================================================================//
//                             kiss_fft_log.h                               //
//==========================================================================//
// [kiss_fft_log.h]
///
#ifndef kiss_fft_log_h
#define kiss_fft_log_h

#define ERROR 1
#define WARNING 2
#define INFO 3
#define DEBUG 4

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#if defined(NDEBUG)
# define KISS_FFT_LOG_MSG(severity, ...) ((void)0)
#else
# define KISS_FFT_LOG_MSG(severity, ...) \
fprintf(stderr, "[" #severity "] " __FILE__ ":" TOSTRING(__LINE__) " "); \
fprintf(stderr, __VA_ARGS__); \
fprintf(stderr, "\n")
#endif

#define KISS_FFT_ERROR(...) KISS_FFT_LOG_MSG(ERROR, __VA_ARGS__)
#define KISS_FFT_WARNING(...) KISS_FFT_LOG_MSG(WARNING, __VA_ARGS__)
#define KISS_FFT_INFO(...) KISS_FFT_LOG_MSG(INFO, __VA_ARGS__)
#define KISS_FFT_DEBUG(...) KISS_FFT_LOG_MSG(DEBUG, __VA_ARGS__)

#endif



//==========================================================================//
//                               kiss_fft.h                                 //
//==========================================================================//
// [kiss_fft.h]
///
#ifndef KISS_FFT_H
#define KISS_FFT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Define KISS_FFT_SHARED macro to properly export symbols
#ifdef KISS_FFT_SHARED
# ifdef _WIN32
#  ifdef KISS_FFT_BUILD
#   define KISS_FFT_API __declspec(dllexport)
#  else
#   define KISS_FFT_API __declspec(dllimport)
#  endif
# else
#  define KISS_FFT_API __attribute__ ((visibility ("default")))
# endif
#else
# define KISS_FFT_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 ATTENTION!
 If you would like a :
 -- a utility that will handle the caching of fft objects
 -- real-only (no imaginary time component ) FFT
 -- a multi-dimensional FFT
 -- a command-line utility to perform ffts
 -- a command-line utility to perform fast-convolution filtering
 
 Then see kfc.h kiss_fftr.h kiss_fftnd.h fftutil.c kiss_fastfir.c
 in the tools/ directory.
 */

/* User may override KISS_FFT_MALLOC and/or KISS_FFT_FREE. */
#ifdef USE_SIMD
# include <xmmintrin.h>
# define kiss_fft_scalar __m128
# ifndef KISS_FFT_MALLOC
#  define KISS_FFT_MALLOC(nbytes) _mm_malloc(nbytes,16)
#  define KISS_FFT_ALIGN_CHECK(ptr)
#  define KISS_FFT_ALIGN_SIZE_UP(size) ((size + 15UL) & ~0xFUL)
# endif
# ifndef KISS_FFT_FREE
#  define KISS_FFT_FREE _mm_free
# endif
#else
# define KISS_FFT_ALIGN_CHECK(ptr)
# define KISS_FFT_ALIGN_SIZE_UP(size) (size)
# ifndef KISS_FFT_MALLOC
#  define KISS_FFT_MALLOC malloc
# endif
# ifndef KISS_FFT_FREE
#  define KISS_FFT_FREE free
# endif
#endif


#ifdef FIXED_POINT
#include <stdint.h>
# if (FIXED_POINT == 32)
#  define kiss_fft_scalar int32_t
# else
#  define kiss_fft_scalar int16_t
# endif
#else
# ifndef kiss_fft_scalar
/*  default is float */
#   define kiss_fft_scalar float
# endif
#endif

typedef struct {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
}kiss_fft_cpx;

typedef struct kiss_fft_state* kiss_fft_cfg;

/*
 *  kiss_fft_alloc
 *
 *  Initialize a FFT (or IFFT) algorithm's cfg/state buffer.
 *
 *  typical usage:      kiss_fft_cfg mycfg=kiss_fft_alloc(1024,0,NULL,NULL);
 *
 *  The return value from fft_alloc is a cfg buffer used internally
 *  by the fft routine or NULL.
 *
 *  If lenmem is NULL, then kiss_fft_alloc will allocate a cfg buffer using malloc.
 *  The returned value should be free()d when done to avoid memory leaks.
 *
 *  The state can be placed in a user supplied buffer 'mem':
 *  If lenmem is not NULL and mem is not NULL and *lenmem is large enough,
 *      then the function places the cfg in mem and the size used in *lenmem
 *      and returns mem.
 *
 *  If lenmem is not NULL and ( mem is NULL or *lenmem is not large enough),
 *      then the function returns NULL and places the minimum cfg
 *      buffer size in *lenmem.
 * */

kiss_fft_cfg KISS_FFT_API kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem);

/*
 * kiss_fft(cfg,in_out_buf)
 *
 * Perform an FFT on a complex input buffer.
 * for a forward FFT,
 * fin should be  f[0] , f[1] , ... ,f[nfft-1]
 * fout will be   F[0] , F[1] , ... ,F[nfft-1]
 * Note that each element is complex and can be accessed like
 f[k].r and f[k].i
 * */
void KISS_FFT_API kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);

/*
 A more generic version of the above function. It reads its input from every Nth sample.
 * */
void KISS_FFT_API kiss_fft_stride(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int fin_stride);

/* If kiss_fft_alloc allocated a buffer, it is one contiguous
 buffer and can be simply free()d when no longer needed*/
#define kiss_fft_free KISS_FFT_FREE

/*
 Cleans up some memory that gets managed internally. Not necessary to call, but it might clean up
 your compiler output to call this before you exit.
 */
void KISS_FFT_API kiss_fft_cleanup(void);


/*
 * Returns the smallest integer k, such that k>=n and k has only "fast" factors (2,3,5)
 */
int KISS_FFT_API kiss_fft_next_fast_size(int n);

/* for real ffts, we need an even size */
#define kiss_fftr_next_fast_size_real(n) \
(kiss_fft_next_fast_size( ((n)+1)>>1)<<1)

#ifdef __cplusplus
}
#endif

#endif



//==========================================================================//
//                              kiss_fftr.h                                 //
//==========================================================================//
// [kiss_fftr.h]
///
#ifndef KISS_FTR_H
#define KISS_FTR_H

#ifdef __cplusplus
extern "C" {
#endif


/*
 
 Real optimized version can save about 45% cpu time vs. complex fft of a real seq.
 
 
 
 */

typedef struct kiss_fftr_state *kiss_fftr_cfg;


kiss_fftr_cfg kiss_fftr_alloc(int nfft,int inverse_fft,void * mem, size_t * lenmem);
/*
 nfft must be even
 
 If you don't care to allocate space, use mem = lenmem = NULL
 */


void kiss_fftr(kiss_fftr_cfg cfg,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata);
/*
 input timedata has nfft scalar points
 output freqdata has nfft/2+1 complex points
 */

void kiss_fftri(kiss_fftr_cfg cfg,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata);
/*
 input freqdata has  nfft/2+1 complex points
 output timedata has nfft scalar points
 */

#define kiss_fftr_free KISS_FFT_FREE

#ifdef __cplusplus
}
#endif
#endif



//==========================================================================//
//                            _kiss_fft_guts.h                              //
//==========================================================================//
// [kiss_fft_guts.h]
///
/* kiss_fft.h
 defines kiss_fft_scalar as either short or a float type
 and defines
 typedef struct { kiss_fft_scalar r; kiss_fft_scalar i; }kiss_fft_cpx; */

#ifndef _kiss_fft_guts_h
#define _kiss_fft_guts_h

//#include "kiss_fft.h"
//#include "kiss_fft_log.h"
#include <limits.h>

#define MAXFACTORS 32
/* e.g. an fft of length 128 has 4 factors
 as far as kissfft is concerned
 4*4*4*2
 */

struct kiss_fft_state{
    int nfft;
    int inverse;
    int factors[2*MAXFACTORS];
    kiss_fft_cpx twiddles[1];
};

/*
 Explanation of macros dealing with complex math:
 
 C_MUL(m,a,b)         : m = a*b
 C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
 C_SUB( res, a,b)     : res = a - b
 C_SUBFROM( res , a)  : res -= a
 C_ADDTO( res , a)    : res += a
 * */
#ifdef FIXED_POINT
#include <stdint.h>
#if (FIXED_POINT==32)
# define FRACBITS 31
# define SAMPPROD int64_t
#define SAMP_MAX INT32_MAX
#define SAMP_MIN INT32_MIN
#else
# define FRACBITS 15
# define SAMPPROD int32_t
#define SAMP_MAX INT16_MAX
#define SAMP_MIN INT16_MIN
#endif

#if defined(CHECK_OVERFLOW)
#  define CHECK_OVERFLOW_OP(a,op,b)  \
if ( (SAMPPROD)(a) op (SAMPPROD)(b) > SAMP_MAX || (SAMPPROD)(a) op (SAMPPROD)(b) < SAMP_MIN ) { \
KISS_FFT_WARNING("overflow (%d " #op" %d) = %ld", (a),(b),(SAMPPROD)(a) op (SAMPPROD)(b)); }
#endif


#   define smul(a,b) ( (SAMPPROD)(a)*(b) )
#   define sround( x )  (kiss_fft_scalar)( ( (x) + (1<<(FRACBITS-1)) ) >> FRACBITS )

#   define S_MUL(a,b) sround( smul(a,b) )

#   define C_MUL(m,a,b) \
do{ (m).r = sround( smul((a).r,(b).r) - smul((a).i,(b).i) ); \
(m).i = sround( smul((a).r,(b).i) + smul((a).i,(b).r) ); }while(0)

#   define DIVSCALAR(x,k) \
(x) = sround( smul(  x, SAMP_MAX/k ) )

#   define C_FIXDIV(c,div) \
do {    DIVSCALAR( (c).r , div);  \
DIVSCALAR( (c).i  , div); }while (0)

#   define C_MULBYSCALAR( c, s ) \
do{ (c).r =  sround( smul( (c).r , s ) ) ;\
(c).i =  sround( smul( (c).i , s ) ) ; }while(0)

#else  /* not FIXED_POINT*/

#   define S_MUL(a,b) ( (a)*(b) )
#define C_MUL(m,a,b) \
do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
(m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
#   define C_FIXDIV(c,div) /* NOOP */
#   define C_MULBYSCALAR( c, s ) \
do{ (c).r *= (s);\
(c).i *= (s); }while(0)
#endif

#ifndef CHECK_OVERFLOW_OP
#  define CHECK_OVERFLOW_OP(a,op,b) /* noop */
#endif

#define  C_ADD( res, a,b)\
do { \
CHECK_OVERFLOW_OP((a).r,+,(b).r)\
CHECK_OVERFLOW_OP((a).i,+,(b).i)\
(res).r=(a).r+(b).r;  (res).i=(a).i+(b).i; \
}while(0)
#define  C_SUB( res, a,b)\
do { \
CHECK_OVERFLOW_OP((a).r,-,(b).r)\
CHECK_OVERFLOW_OP((a).i,-,(b).i)\
(res).r=(a).r-(b).r;  (res).i=(a).i-(b).i; \
}while(0)
#define C_ADDTO( res , a)\
do { \
CHECK_OVERFLOW_OP((res).r,+,(a).r)\
CHECK_OVERFLOW_OP((res).i,+,(a).i)\
(res).r += (a).r;  (res).i += (a).i;\
}while(0)

#define C_SUBFROM( res , a)\
do {\
CHECK_OVERFLOW_OP((res).r,-,(a).r)\
CHECK_OVERFLOW_OP((res).i,-,(a).i)\
(res).r -= (a).r;  (res).i -= (a).i; \
}while(0)


#ifdef FIXED_POINT
#  define KISS_FFT_COS(phase)  floor(.5+SAMP_MAX * cos (phase))
#  define KISS_FFT_SIN(phase)  floor(.5+SAMP_MAX * sin (phase))
#  define HALF_OF(x) ((x)>>1)
#elif defined(USE_SIMD)
#  define KISS_FFT_COS(phase) _mm_set1_ps( cos(phase) )
#  define KISS_FFT_SIN(phase) _mm_set1_ps( sin(phase) )
#  define HALF_OF(x) ((x)*_mm_set1_ps(.5))
#else
#  define KISS_FFT_COS(phase) (kiss_fft_scalar) cos(phase)
#  define KISS_FFT_SIN(phase) (kiss_fft_scalar) sin(phase)
#  define HALF_OF(x) ((x)*((kiss_fft_scalar).5))
#endif

#define  kf_cexp(x,phase) \
do{ \
(x)->r = KISS_FFT_COS(phase);\
(x)->i = KISS_FFT_SIN(phase);\
}while(0)


/* a debugging function */
#define pcpx(c)\
KISS_FFT_DEBUG("%g + %gi\n",(double)((c)->r),(double)((c)->i))


#ifdef KISS_FFT_USE_ALLOCA
// define this to allow use of alloca instead of malloc for temporary buffers
// Temporary buffers are used in two case:
// 1. FFT sizes that have "bad" factors. i.e. not 2,3 and 5
// 2. "in-place" FFTs.  Notice the quotes, since kissfft does not really do an in-place transform.
#include <alloca.h>
#define  KISS_FFT_TMP_ALLOC(nbytes) alloca(nbytes)
#define  KISS_FFT_TMP_FREE(ptr)
#else
#define  KISS_FFT_TMP_ALLOC(nbytes) KISS_FFT_MALLOC(nbytes)
#define  KISS_FFT_TMP_FREE(ptr) KISS_FFT_FREE(ptr)
#endif

#endif



//==========================================================================//
//                              kiss_fftr.c                                 //
//==========================================================================//
// [kiss_fftr.c]
///
struct kiss_fftr_state{
    kiss_fft_cfg substate;
    kiss_fft_cpx * tmpbuf;
    kiss_fft_cpx * super_twiddles;
#ifdef USE_SIMD
    void * pad;
#endif
};

kiss_fftr_cfg kiss_fftr_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem)
{
    KISS_FFT_ALIGN_CHECK(mem)
    
    int i;
    kiss_fftr_cfg st = NULL;
    size_t subsize = 0, memneeded;
    
    if (nfft & 1) {
        KISS_FFT_ERROR("Real FFT optimization must be even.");
        return NULL;
    }
    nfft >>= 1;
    
    kiss_fft_alloc (nfft, inverse_fft, NULL, &subsize);
    memneeded = sizeof(struct kiss_fftr_state) + subsize + sizeof(kiss_fft_cpx) * ( nfft * 3 / 2);
    
    if (lenmem == NULL) {
        st = (kiss_fftr_cfg) KISS_FFT_MALLOC (memneeded);
    } else {
        if (*lenmem >= memneeded)
            st = (kiss_fftr_cfg) mem;
        *lenmem = memneeded;
    }
    if (!st)
        return NULL;
    
    st->substate = (kiss_fft_cfg) (st + 1); /*just beyond kiss_fftr_state struct */
    st->tmpbuf = (kiss_fft_cpx *) (((char *) st->substate) + subsize);
    st->super_twiddles = st->tmpbuf + nfft;
    kiss_fft_alloc(nfft, inverse_fft, st->substate, &subsize);
    
    for (i = 0; i < nfft/2; ++i) {
        double phase =
        -3.14159265358979323846264338327 * ((double) (i+1) / nfft + .5);
        if (inverse_fft)
            phase *= -1;
        kf_cexp (st->super_twiddles+i,phase);
    }
    return st;
}

void kiss_fftr(kiss_fftr_cfg st,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata)
{
    /* input buffer timedata is stored row-wise */
    int k,ncfft;
    kiss_fft_cpx fpnk,fpk,f1k,f2k,tw,tdc;
    
    if ( st->substate->inverse) {
        KISS_FFT_ERROR("kiss fft usage error: improper alloc");
        return;/* The caller did not call the correct function */
    }
    
    ncfft = st->substate->nfft;
    
    /*perform the parallel fft of two real signals packed in real,imag*/
    kiss_fft( st->substate , (const kiss_fft_cpx*)timedata, st->tmpbuf );
    /* The real part of the DC element of the frequency spectrum in st->tmpbuf
     * contains the sum of the even-numbered elements of the input time sequence
     * The imag part is the sum of the odd-numbered elements
     *
     * The sum of tdc.r and tdc.i is the sum of the input time sequence.
     *      yielding DC of input time sequence
     * The difference of tdc.r - tdc.i is the sum of the input (dot product) [1,-1,1,-1...
     *      yielding Nyquist bin of input time sequence
     */
    
    tdc.r = st->tmpbuf[0].r;
    tdc.i = st->tmpbuf[0].i;
    C_FIXDIV(tdc,2);
    CHECK_OVERFLOW_OP(tdc.r ,+, tdc.i);
    CHECK_OVERFLOW_OP(tdc.r ,-, tdc.i);
    freqdata[0].r = tdc.r + tdc.i;
    freqdata[ncfft].r = tdc.r - tdc.i;
#ifdef USE_SIMD
    freqdata[ncfft].i = freqdata[0].i = _mm_set1_ps(0);
#else
    freqdata[ncfft].i = freqdata[0].i = 0;
#endif
    
    for ( k=1;k <= ncfft/2 ; ++k ) {
        fpk    = st->tmpbuf[k];
        fpnk.r =   st->tmpbuf[ncfft-k].r;
        fpnk.i = - st->tmpbuf[ncfft-k].i;
        C_FIXDIV(fpk,2);
        C_FIXDIV(fpnk,2);
        
        C_ADD( f1k, fpk , fpnk );
        C_SUB( f2k, fpk , fpnk );
        C_MUL( tw , f2k , st->super_twiddles[k-1]);
        
        freqdata[k].r = HALF_OF(f1k.r + tw.r);
        freqdata[k].i = HALF_OF(f1k.i + tw.i);
        freqdata[ncfft-k].r = HALF_OF(f1k.r - tw.r);
        freqdata[ncfft-k].i = HALF_OF(tw.i - f1k.i);
    }
}

void kiss_fftri(kiss_fftr_cfg st,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata)
{
    /* input buffer timedata is stored row-wise */
    int k, ncfft;
    
    if (st->substate->inverse == 0) {
        KISS_FFT_ERROR("kiss fft usage error: improper alloc");
        return;/* The caller did not call the correct function */
    }
    
    ncfft = st->substate->nfft;
    
    st->tmpbuf[0].r = freqdata[0].r + freqdata[ncfft].r;
    st->tmpbuf[0].i = freqdata[0].r - freqdata[ncfft].r;
    C_FIXDIV(st->tmpbuf[0],2);
    
    for (k = 1; k <= ncfft / 2; ++k) {
        kiss_fft_cpx fk, fnkc, fek, fok, tmp;
        fk = freqdata[k];
        fnkc.r = freqdata[ncfft - k].r;
        fnkc.i = -freqdata[ncfft - k].i;
        C_FIXDIV( fk , 2 );
        C_FIXDIV( fnkc , 2 );
        
        C_ADD (fek, fk, fnkc);
        C_SUB (tmp, fk, fnkc);
        C_MUL (fok, tmp, st->super_twiddles[k-1]);
        C_ADD (st->tmpbuf[k],     fek, fok);
        C_SUB (st->tmpbuf[ncfft - k], fek, fok);
#ifdef USE_SIMD
        st->tmpbuf[ncfft - k].i *= _mm_set1_ps(-1.0);
#else
        st->tmpbuf[ncfft - k].i *= -1;
#endif
    }
    kiss_fft (st->substate, st->tmpbuf, (kiss_fft_cpx *) timedata);
}



//==========================================================================//
//                                kiss_fft.c                                //
//==========================================================================//
// [kiss_fft.c]
///
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */
static void kf_bfly2(
                     kiss_fft_cpx * Fout,
                     const size_t fstride,
                     const kiss_fft_cfg st,
                     int m
                     )
{
    kiss_fft_cpx * Fout2;
    kiss_fft_cpx * tw1 = st->twiddles;
    kiss_fft_cpx t;
    Fout2 = Fout + m;
    do{
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);
        
        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
        ++Fout2;
        ++Fout;
    }while (--m);
}

static void kf_bfly4(
                     kiss_fft_cpx * Fout,
                     const size_t fstride,
                     const kiss_fft_cfg st,
                     const size_t m
                     )
{
    kiss_fft_cpx *tw1,*tw2,*tw3;
    kiss_fft_cpx scratch[6];
    size_t k=m;
    const size_t m2=2*m;
    const size_t m3=3*m;
    
    
    tw3 = tw2 = tw1 = st->twiddles;
    
    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);
        
        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );
        
        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );
        
        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        ++Fout;
    }while(--k);
}

static void kf_bfly3(
                     kiss_fft_cpx * Fout,
                     const size_t fstride,
                     const kiss_fft_cfg st,
                     size_t m
                     )
{
    size_t k=m;
    const size_t m2 = 2*m;
    kiss_fft_cpx *tw1,*tw2;
    kiss_fft_cpx scratch[5];
    kiss_fft_cpx epi3;
    epi3 = st->twiddles[fstride*m];
    
    tw1=tw2=st->twiddles;
    
    do{
        C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);
        
        C_MUL(scratch[1],Fout[m] , *tw1);
        C_MUL(scratch[2],Fout[m2] , *tw2);
        
        C_ADD(scratch[3],scratch[1],scratch[2]);
        C_SUB(scratch[0],scratch[1],scratch[2]);
        tw1 += fstride;
        tw2 += fstride*2;
        
        Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
        Fout[m].i = Fout->i - HALF_OF(scratch[3].i);
        
        C_MULBYSCALAR( scratch[0] , epi3.i );
        
        C_ADDTO(*Fout,scratch[3]);
        
        Fout[m2].r = Fout[m].r + scratch[0].i;
        Fout[m2].i = Fout[m].i - scratch[0].r;
        
        Fout[m].r -= scratch[0].i;
        Fout[m].i += scratch[0].r;
        
        ++Fout;
    }while(--k);
}

static void kf_bfly5(
                     kiss_fft_cpx * Fout,
                     const size_t fstride,
                     const kiss_fft_cfg st,
                     int m
                     )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];
    
    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;
    
    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;
        
        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);
        
        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);
        
        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;
        
        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);
        
        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);
        
        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);
        
        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);
        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);
        
        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);
        
        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(
                            kiss_fft_cpx * Fout,
                            const size_t fstride,
                            const kiss_fft_cfg st,
                            int m,
                            int p
                            )
{
    int u,k,q1,q;
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx t;
    int Norig = st->nfft;
    
    kiss_fft_cpx * scratch = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC(sizeof(kiss_fft_cpx)*p);
    if (scratch == NULL){
        KISS_FFT_ERROR("Memory allocation failed.");
        return;
    }
    
    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
            C_FIXDIV(scratch[q1],p);
            k += m;
        }
        
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratch[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
    KISS_FFT_TMP_FREE(scratch);
}

static
void kf_work(
             kiss_fft_cpx * Fout,
             const kiss_fft_cpx * f,
             const size_t fstride,
             int in_stride,
             int * factors,
             const kiss_fft_cfg st
             )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const kiss_fft_cpx * Fout_end = Fout + p*m;
    
#ifdef _OPENMP
    // use openmp extensions at the
    // top-level (not recursive)
    if (fstride==1 && p<=5 && m!=1)
    {
        int k;
        
        // execute the p different work units in different threads
#       pragma omp parallel for
        for (k=0;k<p;++k)
        kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        // all threads have joined by this point
        
        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break;
            case 4: kf_bfly4(Fout,fstride,st,m); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break;
            default: kf_bfly_generic(Fout,fstride,st,m,p); break;
        }
        return;
    }
#endif
    
    if (m==1) {
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }
    
    Fout=Fout_beg;
    
    // recombine the p smaller DFTs
    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break;
        case 4: kf_bfly4(Fout,fstride,st,m); break;
        case 5: kf_bfly5(Fout,fstride,st,m); break;
        default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
 where
 p[i] * m[i] = m[i-1]
 m0 = n                  */
static
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );
    
    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
kiss_fft_cfg kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem )
{
    KISS_FFT_ALIGN_CHECK(mem)
    
    kiss_fft_cfg st=NULL;
    size_t memneeded = KISS_FFT_ALIGN_SIZE_UP(sizeof(struct kiss_fft_state)
                                              + sizeof(kiss_fft_cpx)*(nfft-1)); /* twiddle factors*/
    
    if ( lenmem==NULL ) {
        st = ( kiss_fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;
        
        for (i=0;i<nfft;++i) {
            const double pi=3.141592653589793238462643383279502884197169399375105820974944;
            double phase = -2*pi*i / nfft;
            if (st->inverse)
                phase *= -1;
            kf_cexp(st->twiddles+i, phase );
        }
        
        kf_factor(nfft,st->factors);
    }
    return st;
}


void kiss_fft_stride(kiss_fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        //NOTE: this is not really an in-place FFT algorithm.
        //It just performs an out-of-place FFT into a temp buffer
        if (fout == NULL){
            KISS_FFT_ERROR("fout buffer NULL.");
            return;
        }
        
        kiss_fft_cpx * tmpbuf = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC( sizeof(kiss_fft_cpx)*st->nfft);
        if (tmpbuf == NULL){
            KISS_FFT_ERROR("Memory allocation error.");
            return;
        }
        
        
        
        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
        KISS_FFT_TMP_FREE(tmpbuf);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    kiss_fft_stride(cfg,fin,fout,1);
}


void kiss_fft_cleanup(void)
{
    // nothing needed any more
}

int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}



//==========================================================================//
//                               CASTOR TOOLS                               //
//==========================================================================//
namespace castor
{

//==========================================================================
// [ckiss]
/// Complex class for kissfft interface.
class ckiss : public kiss_fft_cpx
{
public:
    ckiss() {r=0; i=0;};
    template<typename S>
    ckiss(S v, S w=0) {r=v; i=w;};
    template<typename S>
    ckiss(std::complex<S> v) {r=v.real(); i=v.imag();};
    template<typename S>
    operator std::complex<S>() const {return std::complex<S>(r,i);};
    ckiss& operator*=(ckiss const& c)
    {
        float t = r;
        r=r*c.r - i*c.i;
        i=t*c.i + i*c.r;
        return (*this);
    }
};
std::ostream& operator<<(std::ostream& flux, ckiss const& c)
{
    flux << "(" << c.r << "," << c.i << ")"; return flux;
}
inline void disp(matrix<ckiss>const& A, int info=2, std::ostream& flux=std::cout, std::size_t m=3, std::size_t n=3)
{
    disp<std::complex<float>>(A,info,flux,m,n);
}

//==========================================================================
// [fft]
/// Fast Fourier Transform of an array.
///
/// Y = fft(X) computes the discrete Fourier transform (DFT) of X using a
/// fast Fourier transform (FFT) algorithm.
///
/// Y = fft(X,n) returns the n-point DFT. If no value is specified, Y is the
/// fft of all values of X, in linear indexing.
///
/// Y = fft(X,n,dim) returns the Fourier transform along the dimension dim.
/// For example, if X is a matrix, then fft(X,n,2) returns the n-point Fourier
/// transform of each row.
///
/// Y = fft(X,n,dim,flag) specify is inverse fft is computed.
///
/// WARNING : This is a single-precision floating point implementation.
///
/// \code{.cpp}
///    matrix<float> X = eye(1,4);
///    auto Y = fft(X);
///    disp(Y);
/// \endcode
///
// \see dft, ifft.
matrix<std::complex<float>> fft(matrix<ckiss> Xc, std::size_t n=0, int dim=0, bool isinverse=false)
{
    if ((dim<0) || (dim>2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimension argument must be 0 for all, 1 for rows or 2 for columns.");
    }
    if (n==0) {n = Xc.size(dim);}
    matrix<ckiss> Yc;
    kiss_fft_cfg cfg = kiss_fft_alloc((int)n, isinverse, NULL, NULL);
    if (dim==0)
    {
        Xc.reshape(1,numel(Xc));
        Xc.resize(1,n,0);
        Yc = zeros<ckiss>(size(Xc));
        kiss_fft(cfg,&Xc(0),&Yc(0));
    }
    else if (dim==1)
    {
        Xc.resize(n,size(Xc,2),0);
        Xc = transpose(Xc);
        Yc = zeros<ckiss>(size(Xc));
        for (std::size_t j=0; j<size(Xc,1); ++j)
        {
            kiss_fft(cfg,&Xc(j,0),&Yc(j,0));
        }
        Yc = transpose(Yc);
    }
    else if (dim==2)
    {
        Xc.resize(size(Xc,1),n,0);
        Yc = zeros<ckiss>(size(Xc));
        for (std::size_t i=0; i<size(Xc,1); ++i)
        {
            kiss_fft(cfg,&Xc(i,0),&Yc(i,0));
        }
    }
    kiss_fft_free(cfg);
    if (isinverse)
    {
        for (std::size_t l=0; l<numel(Yc); ++l)
        {
            Yc(l).r /= n;
            Yc(l).i /= n;
        }
    }
    return Yc;
}

//==========================================================================
// [fftconv]
/// Convolution product and polynomial multiplication using FFT.
///
/// C = fftconv(A,B) convolves array A and B using Fast Fourier Transform
/// that guarantee efficient execution. The resulting array is vector
/// of length NUMEL(A)+NUMEL(B)-1. If A and B are vectors of polynomial
/// coefficients, convolving them is equivalent to multiplying the two
/// polynomials.
///
/// C = fftconv(A,B,DIM) convolve along dimension DIM if A and B have
/// compatible size.
///
/// \code{.cpp}
///    matrix<float> A = eye(3,4);
///    matrix<float> B = rand(3,10);
///    disp(fftconv(A,B,2);
/// \endcode
///
// \see conv, fft.
template<typename S>
auto fftconv(matrix<S>const& A, matrix<S>const& B, int dim=0)
{
    using T = std::complex<decltype(std::abs(A(0)))>;
    std::size_t n = A.size(dim) + B.size(dim)-1;
    std::size_t nfft = std::pow(2,(std::size_t)std::log2(n)+1);
    auto Af = fft(A,nfft,dim);
    auto Bf = fft(B,nfft,dim);
    Af *= Bf;
    matrix<T> C = ifft(Af,0,dim);
    if (dim==0)
    {
        C = eval(C(range(0,n)));
    }
    else if (dim==1)
    {
        C = eval(C(range(0,n),col(C)));
    }
    else if (dim==2)
    {
        C = eval(C(row(C),range(0,n)));
    }
    return C;
}

//==========================================================================
// [ifft]
/// Inverse Fast Fourier Transform of an array.
///
/// Y = ifft(X) computes the inverse discrete Fourier transform (IDFT) of X
/// using a fast Fourier transform (FFT) algorithm.
///
/// Y = ifft(X,n) returns the n-point IDFT. If no value is specified, Y is the
/// ifft of all values of X, in linear indexing.
///
/// Y = ifft(X,n,dim) returns the inverse Fourier transform along the dimension dim.
/// For example, if X is a matrix, then ifft(X,n,2) returns the n-point inverse
/// Fourier transform of each row.
///
/// WARNING : This is a single-precision floating point implementation.
///
/// \code{.cpp}
///    matrix<float> X = eye(1,512);
///    auto Y = ifft(X);
///    disp(Y);
/// \endcode
///
// \see idft, fft.
inline matrix<std::complex<float>> ifft(matrix<ckiss> Xc, std::size_t n=0, int dim=0)
{
    return fft(Xc, n, dim, true);
}


}
