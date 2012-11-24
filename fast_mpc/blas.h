#ifndef BLAS_H
#define BLAS_H
/** \file   blas.h
 *  \brief  Header file for Fortran BLAS.
 */


#define F77_CALL(x)  x ## _
#define F77_NAME(x)  F77_CALL(x)

#ifdef  __cplusplus
extern "C" {
#endif

/* Level 1 BLAS */

extern double F77_NAME(dasum)(const size_t *n, const double *dx, const size_t *incx);
extern void   F77_NAME(daxpy)(const size_t *n, const double *alpha,
                              const double *dx, const size_t *incx,
                              double *dy, const size_t *incy);
extern void   F77_NAME(dcopy)(const size_t *n, const double *dx, const size_t *incx,
                              double *dy, const size_t *incy);
extern double F77_NAME(ddot) (const size_t *n, const double *dx, const size_t *incx,
                              const double *dy, const size_t *incy);
extern double F77_NAME(dnrm2)(const size_t *n, const double *dx, const size_t *incx);
extern void   F77_NAME(drot) (const size_t *n, double *dx, const size_t *incx,
                              double *dy, const size_t *incy, const double *c,
                              const double *s);
extern void   F77_NAME(drotg)(const double *a, const double *b, double *c, 
                              double *s);
extern void   F77_NAME(drotm)(const size_t *n, double *dx, const size_t *incx,
                              double *dy, const size_t *incy,const double *dparam);
extern void   F77_NAME(drotmg)(const double *dd1, const double *dd2, 
                               const double *dx1, const double *dy1, 
                               double *param);
extern void   F77_NAME(dscal)(const size_t *n, const double *alpha, double *dx,
                              const size_t *incx);
extern void   F77_NAME(dswap)(const size_t *n, double *dx, const size_t *incx,
                              double *dy, const size_t *incy);
extern size_t    F77_NAME(idamax)(const size_t *n, const double *dx, const size_t *incx);

/* Level 2 BLAS */

extern void   F77_NAME(dgbmv)(const char *trans, const size_t *m, const size_t *n,
                              const size_t *kl,const size_t *ku, const double *alpha,
                              const double *a, const size_t *lda, const double *x,
                              const size_t *incx, const double *beta, double *y,
                              const size_t *incy);
extern void   F77_NAME(dgemv)(const char *trans, const size_t *m, const size_t *n,
                              const double *alpha, const double *a,
                              const size_t *lda, const double *x, const size_t *incx,
                              const double *beta, double *y, const size_t *incy);
extern void   F77_NAME(dsbmv)(const char *uplo, const size_t *n, const size_t *k,
                              const double *alpha, const double *a,
                              const size_t *lda, const double *x, const size_t *incx,
                              const double *beta, double *y, const size_t *incy);
extern void   F77_NAME(dspmv)(const char *uplo, const size_t *n,
                              const double *alpha, const double *ap,
                              const double *x, const size_t *incx,
                              const double *beta, double *y, const size_t *incy);
extern void   F77_NAME(dsymv)(const char *uplo, const size_t *n,
                              const double *alpha, const double *a,
                              const size_t *lda, const double *x, const size_t *incx,
                              const double *beta, double *y, const size_t *incy);
extern void   F77_NAME(dtbmv)(const char *uplo, const char *trans,
                              const char *diag, const size_t *n, const size_t *k,
                              const double *a, const size_t *lda,
                              double *x, const size_t *incx);
extern void   F77_NAME(dtpmv)(const char *uplo, const char *trans,
                              const char *diag, const size_t *n, const double *ap,
                              double *x, const size_t *incx);
extern void   F77_NAME(dtrmv)(const char *uplo, const char *trans,
                              const char *diag, const size_t *n, const double *a,
                              const size_t *lda, double *x, const size_t *incx);
extern void   F77_NAME(dtbsv)(const char *uplo, const char *trans,
                              const char *diag, const size_t *n, const size_t *k,
                              const double *a, const size_t *lda,
                              double *x, const size_t *incx);
extern void   F77_NAME(dtpsv)(const char *uplo, const char *trans,
                              const char *diag, const size_t *n,
                              const double *ap, double *x, const size_t *incx);
extern void   F77_NAME(dtrsv)(const char *uplo, const char *trans,
                              const char *diag, const size_t *n,
                              const double *a, const size_t *lda,
                              double *x, const size_t *incx);
extern void   F77_NAME(dger) (const size_t *m, const size_t *n, const double *alpha,
                              double *x, const size_t *incx,
                              double *y, const size_t *incy,
                              double *a, const size_t *lda);
extern void   F77_NAME(dsyr) (const char *uplo, const size_t *n,
                              const double *alpha, const double *x,
                              const size_t *incx, double *a, const size_t *lda);
extern void   F77_NAME(dspr) (const char *uplo, const size_t *n,
                              const double *alpha, const double *x,
                              const size_t *incx, double *ap);
extern void   F77_NAME(dsyr2)(const char *uplo, const size_t *n, 
                              const double *alpha, const double *x,
                              const size_t *incx, const double *y, const size_t *incy,
                              double *a, const size_t *lda);
extern void   F77_NAME(dspr2)(const char *uplo, const size_t *n,
                              const double *alpha, const double *x,
                              const size_t *incx, const double *y,
                              const size_t *incy, double *ap);

/* Level 3 BLAS */

extern void   F77_NAME(dgemm)(const char *transa, const char *transb,
                              const size_t *m, const size_t *n, const size_t *k,
                              const double *alpha, const double *a,
                              const size_t *lda, const double *b, const size_t *ldb,
                              const double *beta, double *c, const size_t *ldc);
extern void   F77_NAME(dtrsm)(const char *side, const char *uplo,
                              const char *transa, const char *diag,
                              const size_t *m, const size_t *n, const double *alpha,
                              const double *a, const size_t *lda,
                              double *b, const size_t *ldb);
extern void   F77_NAME(dtrmm)(const char *side, const char *uplo,
                              const char *transa, const char *diag,
                              const size_t *m, const size_t *n, const double *alpha,
                              const double *a, const size_t *lda,
                              double *b, const size_t *ldb);
extern void   F77_NAME(dsymm)(const char *side, const char *uplo, const size_t *m,
                              const size_t *n, const double *alpha,
                              const double *a, const size_t *lda,
                              const double *b, const size_t *ldb,
                              const double *beta, double *c, const size_t *ldc);
extern void   F77_NAME(dsyrk)(const char *uplo, const char *trans,
                              const size_t *n, const size_t *k,
                              const double *alpha, const double *a,
                              const size_t *lda, const double *beta,
                              double *c, const size_t *ldc);
extern void   F77_NAME(dsyr2k)(const char *uplo, const char *trans,
                               const size_t *n, const size_t *k,
                               const double *alpha, const double *a,
                               const size_t *lda, const double *b, const size_t *ldb,
                               const double *beta, double *c, const size_t *ldc);

#ifdef  __cplusplus
}
#endif

#endif /* BLAS_H */
