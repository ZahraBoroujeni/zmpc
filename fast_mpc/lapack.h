#ifndef LAPACK_H
#define LAPACK_H

#include "blas.h"

#ifdef  __cplusplus
extern "C" {
#endif

extern void F77_NAME(dbdsqr)(const char* uplo, const size_t* n, const size_t* ncvt,
		 const size_t* nru, const size_t* ncc, double* d, double* e,
		 double* vt, const size_t* ldvt, double* u, const size_t* ldu,
		 double* c, const size_t* ldc, double* work, size_t* info);
extern void F77_NAME(ddisna)(const char* job, const size_t* m, const size_t* n,
		 double* d, double* sep, size_t* info);

extern void F77_NAME(dgbbrd)(const char* vect, const size_t* m, const size_t* n,
		 const size_t* ncc, const size_t* kl, const size_t* ku,
		 double* ab, const size_t* ldab,
		 double* d, double* e, double* q,
		 const size_t* ldq, double* pt, const size_t* ldpt,
		 double* c, const size_t* ldc,
		 double* work, size_t* info);
extern void F77_NAME(dgbcon)(const char* norm, const size_t* n, const size_t* kl,
		 const size_t* ku, double* ab, const size_t* ldab,
		 size_t* ipiv, const double* anorm, double* rcond,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgbequ)(const size_t* m, const size_t* n, const size_t* kl, const size_t* ku,
		 double* ab, const size_t* ldab, double* r, double* c,
		 double* rowcnd, double* colcnd, double* amax, size_t* info);
extern void F77_NAME(dgbrfs)(const char* trans, const size_t* n, const size_t* kl,
		 const size_t* ku, const size_t* nrhs, double* ab,
		 const size_t* ldab, double* afb, const size_t* ldafb,
		 size_t* ipiv, double* b, const size_t* ldb,
		 double* x, const size_t* ldx, double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgbsv)(const size_t* n, const size_t* kl,const size_t* ku,
		const size_t* nrhs, double* ab, const size_t* ldab,
		size_t* ipiv, double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dgbsvx)(const size_t* fact, const char* trans,
		 const size_t* n, const size_t* kl,const size_t* ku,
		 const size_t* nrhs, double* ab, const size_t* ldab,
		 double* afb, const size_t* ldafb, size_t* ipiv,
		 const char* equed, double* r, double* c, 
		 double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgbtf2)(const size_t* m, const size_t* n, const size_t* kl,const size_t* ku,
		 double* ab, const size_t* ldab, size_t* ipiv, size_t* info);
extern void F77_NAME(dgbtrf)(const size_t* m, const size_t* n, const size_t* kl,const size_t* ku,
    		  double* ab, const size_t* ldab, size_t* ipiv, size_t* info);
extern void F77_NAME(dgbtrs)(const char* trans, const size_t* n,
		 const size_t* kl, const size_t* ku, const size_t* nrhs,
		 const double* ab, const size_t* ldab, const size_t* ipiv,
		 double* b, const size_t* ldb, size_t* info);

extern void F77_NAME(dgebak)(const char* job, const char* side, const size_t* n,
		 const size_t* ilo, const size_t* ihi, double* scale,
		 const size_t* m, double* v, const size_t* ldv, size_t* info);
extern void F77_NAME(dgebal)(const char* job, const size_t* n, double* a, const size_t* lda,
    		  size_t* ilo, size_t* ihi, double* scale, size_t* info);
extern void F77_NAME(dgebd2)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* d, double* e, double* tauq, double* taup,
		 double* work, size_t* info);
extern void F77_NAME(dgebrd)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* d, double* e, double* tauq, double* taup,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgecon)(const char* norm, const size_t* n,
		 const double* a, const size_t* lda,
		 const double* anorm, double* rcond,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgeequ)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* r, double* c, double* rowcnd, double* colcnd,
		 double* amax, size_t* info);
extern void F77_NAME(dgees)(const char* jobvs, const char* sort,
		size_t (*select)(const double*, const double*),
		const size_t* n, double* a, const size_t* lda,
		size_t* sdim, double* wr, double* wi,
		double* vs, const size_t* ldvs,
		double* work, const size_t* lwork, size_t* bwork, size_t* info);
extern void F77_NAME(dgeesx)(const char* jobvs, const char* sort,
		 size_t (*select)(const double*, const double*),
		 const char* sense, const size_t* n, double* a,
		 const size_t* lda, size_t* sdim, double* wr, double* wi,
		 double* vs, const size_t* ldvs, double* rconde,
		 double* rcondv, double* work, const size_t* lwork,
		 size_t* iwork, const size_t* liwork, size_t* bwork, size_t* info);
extern void F77_NAME(dgeev)(const char* jobvl, const char* jobvr,
		const size_t* n, double* a, const size_t* lda,
		double* wr, double* wi, double* vl, const size_t* ldvl,
		double* vr, const size_t* ldvr,
		double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgeevx)(const char* balanc, const char* jobvl, const char* jobvr,
		 const char* sense, const size_t* n, double* a, const size_t* lda,
		 double* wr, double* wi, double* vl, const size_t* ldvl,
		 double* vr, const size_t* ldvr, size_t* ilo, size_t* ihi,
		 double* scale, double* abnrm, double* rconde, double* rcondv,
		 double* work, const size_t* lwork, size_t* iwork, size_t* info);
extern void F77_NAME(dgegv)(const char* jobvl, const char* jobvr,
		const size_t* n, double* a, const size_t* lda,
		double* b, const size_t* ldb,
		double* alphar, double* alphai,
		const double* beta, double* vl, const size_t* ldvl,
		double* vr, const size_t* ldvr,
		double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgehd2)(const size_t* n, const size_t* ilo, const size_t* ihi,
		 double* a, const size_t* lda, double* tau,
		 double* work, size_t* info);
extern void F77_NAME(dgehrd)(const size_t* n, const size_t* ilo, const size_t* ihi,
		 double* a, const size_t* lda, double* tau,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgelq2)(const size_t* m, const size_t* n,
		 double* a, const size_t* lda, double* tau,
		 double* work, size_t* info);
extern void F77_NAME(dgelqf)(const size_t* m, const size_t* n,
		 double* a, const size_t* lda, double* tau,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgels)(const char* trans, const size_t* m, const size_t* n,
		const size_t* nrhs, double* a, const size_t* lda,
		double* b, const size_t* ldb,
		double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgelss)(const size_t* m, const size_t* n, const size_t* nrhs,
		 double* a, const size_t* lda, double* b, const size_t* ldb,
		 double* s, double* rcond, size_t* rank,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgelsy)(const size_t* m, const size_t* n, const size_t* nrhs,
		 double* a, const size_t* lda, double* b, const size_t* ldb,
		 size_t* jpvt, const double* rcond, size_t* rank,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgeql2)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* tau, double* work, size_t* info);
extern void F77_NAME(dgeqlf)(const size_t* m, const size_t* n,
		 double* a, const size_t* lda, double* tau,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgeqp3)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 size_t* jpvt, double* tau, double* work, const size_t* lwork,
		 size_t* info);
extern void F77_NAME(dgeqpf)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 size_t* jpvt, double* tau, double* work, size_t* info);
extern void F77_NAME(dgeqr2)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* tau, double* work, size_t* info);
extern void F77_NAME(dgeqrf)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* tau, double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgerfs)(const char* trans, const size_t* n, const size_t* nrhs,
		 double* a, const size_t* lda, double* af, const size_t* ldaf,
		 size_t* ipiv, double* b, const size_t* ldb,
		 double* x, const size_t* ldx, double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgerq2)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* tau, double* work, size_t* info);
extern void F77_NAME(dgerqf)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 double* tau, double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgesv)(const size_t* n, const size_t* nrhs, double* a, const size_t* lda,
		size_t* ipiv, double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dgesvd)(const char* jobu, const char* jobvt, const size_t* m,
		 const size_t* n, double* a, const size_t* lda, double* s,
		 double* u, const size_t* ldu, double* vt, const size_t* ldvt,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgesvx)(const size_t* fact, const char* trans, const size_t* n,
		 const size_t* nrhs, double* a, const size_t* lda,
		 double* af, const size_t* ldaf, size_t* ipiv,
		 char *equed, double* r, double* c,
		 double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgetf2)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 size_t* ipiv, size_t* info);
extern void F77_NAME(dgetrf)(const size_t* m, const size_t* n, double* a, const size_t* lda,
		 size_t* ipiv, size_t* info);
extern void F77_NAME(dgetri)(const size_t* n, double* a, const size_t* lda,
		 size_t* ipiv, double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgetrs)(const char* trans, const size_t* n, const size_t* nrhs,
		 const double* a, const size_t* lda, const size_t* ipiv,
		 double* b, const size_t* ldb, size_t* info);

extern void F77_NAME(dggbak)(const char* job, const char* side,
		 const size_t* n, const size_t* ilo, const size_t* ihi,
		 double* lscale, double* rscale, const size_t* m,
		 double* v, const size_t* ldv, size_t* info);
extern void F77_NAME(dggbal)(const char* job, const size_t* n, double* a, const size_t* lda,
		 double* b, const size_t* ldb, size_t* ilo, size_t* ihi,
		 double* lscale, double* rscale, double* work, size_t* info);
extern void F77_NAME(dgges)(const char* jobvsl, const char* jobvsr, const char* sort,
		size_t (*delztg)(double*, double*, double*),
		const size_t* n, double* a, const size_t* lda,
		double* b, const size_t* ldb, double* alphar,
		double* alphai, const double* beta,
		double* vsl, const size_t* ldvsl,
		double* vsr, const size_t* ldvsr,
		double* work, const size_t* lwork, size_t* bwork, size_t* info);

extern void F77_NAME(dggglm)(const size_t* n, const size_t* m, const size_t* p,
		 double* a, const size_t* lda, double* b, const size_t* ldb,
		 double* d, double* x, double* y,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dgghrd)(const char* compq, const char* compz, const size_t* n,
		 const size_t* ilo, const size_t* ihi, double* a, const size_t* lda,
		 double* b, const size_t* ldb, double* q, const size_t* ldq,
		 double* z, const size_t* ldz, size_t* info);
extern void F77_NAME(dgglse)(const size_t* m, const size_t* n, const size_t* p,
		 double* a, const size_t* lda,
		 double* b, const size_t* ldb,
		 double* c, double* d, double* x,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dggqrf)(const size_t* n, const size_t* m, const size_t* p,
		 double* a, const size_t* lda, double* taua,
		 double* b, const size_t* ldb, double* taub,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dggrqf)(const size_t* m, const size_t* p, const size_t* n,
		 double* a, const size_t* lda, double* taua,
		 double* b, const size_t* ldb, double* taub,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dggsvd)(const char* jobu, const char* jobv, const char* jobq,
		 const size_t* m, const size_t* n, const size_t* p,
		 const size_t* k, const size_t* l,
		 double* a, const size_t* lda,
		 double* b, const size_t* ldb,
		 const double* alpha, const double* beta,
		 double* u, const size_t* ldu,
		 double* v, const size_t* ldv,
		 double* q, const size_t* ldq,
		 double* work, size_t* iwork, size_t* info);

extern void F77_NAME(dgtcon)(const char* norm, const size_t* n, double* dl, double* d,
		 double* du, double* du2, size_t* ipiv, const double* anorm,
		 double* rcond, double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgtrfs)(const char* trans, const size_t* n, const size_t* nrhs,
		 double* dl, double* d, double* du, double* dlf,
		 double* df, double* duf, double* du2,
		 size_t* ipiv, double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgtsv)(const size_t* n, const size_t* nrhs,
		double* dl, double* d, double* du,
		double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dgtsvx)(const size_t* fact, const char* trans,
		 const size_t* n, const size_t* nrhs,
		 double* dl, double* d, double* du,
		 double* dlf, double* df, double* duf,
		 double* du2, size_t* ipiv,
		 double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dgttrf)(const size_t* n, double* dl, double* d,
		 double* du, double* du2, size_t* ipiv, size_t* info);
extern void F77_NAME(dgttrs)(const char* trans, const size_t* n, const size_t* nrhs,
		 double* dl, double* d, double* du, double* du2,
		 size_t* ipiv, double* b, const size_t* ldb, size_t* info);

extern void F77_NAME(dopgtr)(const char* uplo, const size_t* n,
		 const double* ap, const double* tau,
		 double* q, const size_t* ldq,
		 double* work, size_t* info);
extern void F77_NAME(dopmtr)(const char* side, const char* uplo,
		 const char* trans, const size_t* m, const size_t* n,
		 const double* ap, const double* tau,
		 double* c, const size_t* ldc,
		 double* work, size_t* info);
extern void F77_NAME(dorg2l)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda,
		 const double* tau, double* work, size_t* info);
extern void F77_NAME(dorg2r)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda,
		 const double* tau, double* work, size_t* info);
extern void F77_NAME(dorgbr)(const char* vect, const size_t* m,
		 const size_t* n, const size_t* k,
		 double* a, const size_t* lda,
		 const double* tau, double* work,
		 const size_t* lwork, size_t* info);
extern void F77_NAME(dorghr)(const size_t* n, const size_t* ilo, const size_t* ihi,
		 double* a, const size_t* lda, const double* tau,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dorgl2)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda, const double* tau,
		 double* work, size_t* info);
extern void F77_NAME(dorglq)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda,
		 const double* tau, double* work,
		 const size_t* lwork, size_t* info);
extern void F77_NAME(dorgql)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda,
		 const double* tau, double* work,
		 const size_t* lwork, size_t* info);
extern void F77_NAME(dorgqr)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda, const double* tau,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dorgr2)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda, const double* tau,
		 double* work, size_t* info);
extern void F77_NAME(dorgrq)(const size_t* m, const size_t* n, const size_t* k,
		 double* a, const size_t* lda, const double* tau,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dorgtr)(const char* uplo, const size_t* n,
		 double* a, const size_t* lda, const double* tau,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dorm2l)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda,
		 const double* tau, double* c, const size_t* ldc,
		 double* work, size_t* info);
extern void F77_NAME(dorm2r)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda, const double* tau,
		 double* c, const size_t* ldc, double* work, size_t* info);
extern void F77_NAME(dormbr)(const char* vect, const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda, const double* tau,
		 double* c, const size_t* ldc,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dormhr)(const char* side, const char* trans, const size_t* m,
		 const size_t* n, const size_t* ilo, const size_t* ihi,
		 const double* a, const size_t* lda, const double* tau,
		 double* c, const size_t* ldc,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dorml2)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda, const double* tau,
		 double* c, const size_t* ldc, double* work, size_t* info);
extern void F77_NAME(dormlq)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda,
		 const double* tau, double* c, const size_t* ldc,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dormql)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda,
		 const double* tau, double* c, const size_t* ldc,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dormqr)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda,
		 const double* tau, double* c, const size_t* ldc,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dormr2)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda,
		 const double* tau, double* c, const size_t* ldc,
		 double* work, size_t* info);
extern void F77_NAME(dormrq)(const char* side, const char* trans,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* a, const size_t* lda,
		 const double* tau, double* c, const size_t* ldc,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dormtr)(const char* side, const char* uplo,
		 const char* trans, const size_t* m, const size_t* n,
		 const double* a, const size_t* lda,
		 const double* tau, double* c, const size_t* ldc,
		 double* work, const size_t* lwork, size_t* info);

extern void F77_NAME(dpbcon)(const char* uplo, const size_t* n, const size_t* kd,
		 const double* ab, const size_t* ldab,
		 const double* anorm, double* rcond,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dpbequ)(const char* uplo, const size_t* n, const size_t* kd,
		 const double* ab, const size_t* ldab,
		 double* s, double* scond, double* amax, size_t* info);
extern void F77_NAME(dpbrfs)(const char* uplo, const size_t* n,
		 const size_t* kd, const size_t* nrhs,
		 const double* ab, const size_t* ldab,
		 const double* afb, const size_t* ldafb,
		 const double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dpbstf)(const char* uplo, const size_t* n, const size_t* kd,
		 double* ab, const size_t* ldab, size_t* info);
extern void F77_NAME(dpbsv)(const char* uplo, const size_t* n,
		const size_t* kd, const size_t* nrhs,
		double* ab, const size_t* ldab,
		double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dpbsvx)(const size_t* fact, const char* uplo, const size_t* n,
		 const size_t* kd, const size_t* nrhs,
		 double* ab, const size_t* ldab,
		 double* afb, const size_t* ldafb,
		 char* equed, double* s,
		 double* b, const size_t* ldb,
		 double* x, const size_t* ldx, double* rcond,
		 double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dpbtf2)(const char* uplo, const size_t* n, const size_t* kd,
		 double* ab, const size_t* ldab, size_t* info);
extern void F77_NAME(dpbtrf)(const char* uplo, const size_t* n, const size_t* kd,
		 double* ab, const size_t* ldab, size_t* info);
extern void F77_NAME(dpbtrs)(const char* uplo, const size_t* n,
		 const size_t* kd, const size_t* nrhs,
		 const double* ab, const size_t* ldab,
		 double* b, const size_t* ldb, size_t* info);

extern void F77_NAME(dpocon)(const char* uplo, const size_t* n,
		 const double* a, const size_t* lda,
		 const double* anorm, double* rcond,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dpoequ)(const size_t* n, const double* a, const size_t* lda,
		 double* s, double* scond, double* amax, size_t* info);
extern void F77_NAME(dporfs)(const char* uplo, const size_t* n, const size_t* nrhs,
		 const double* a, const size_t* lda,
		 const double* af, const size_t* ldaf,
		 const double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dposv)(const char* uplo, const size_t* n, const size_t* nrhs,
		double* a, const size_t* lda,
		double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dposvx)(const size_t* fact, const char* uplo,
		 const size_t* n, const size_t* nrhs,
		 double* a, const size_t* lda,
		 double* af, const size_t* ldaf, char* equed,
		 double* s, double* b, const size_t* ldb,
		 double* x, const size_t* ldx, double* rcond,
		 double* ferr, double* berr, double* work,
		 size_t* iwork, size_t* info);
extern void F77_NAME(dpotf2)(const char* uplo, const size_t* n,
		 double* a, const size_t* lda, size_t* info);
extern void F77_NAME(dpotrf)(const char* uplo, const size_t* n,
		 double* a, const size_t* lda, size_t* info);
extern void F77_NAME(dpotri)(const char* uplo, const size_t* n,
		 double* a, const size_t* lda, size_t* info);
extern void F77_NAME(dpotrs)(const char* uplo, const size_t* n,
		 const size_t* nrhs, const double* a, const size_t* lda,
		 double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dppcon)(const char* uplo, const size_t* n,
		 const double* ap, const double* anorm, double* rcond,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dppequ)(const char* uplo, const size_t* n,
		 const double* ap, double* s, double* scond,
		 double* amax, size_t* info);

extern void F77_NAME(dpprfs)(const char* uplo, const size_t* n, const size_t* nrhs,
		 const double* ap, const double* afp,
		 const double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dppsv)(const char* uplo, const size_t* n,
		const size_t* nrhs, const double* ap,
		double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dppsvx)(const size_t* fact, const char* uplo,
		 const size_t* n, const size_t* nrhs, double* ap,
		 double* afp, char* equed, double* s,
		 double* b, const size_t* ldb,
		 double* x, const size_t* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dpptrf)(const char* uplo, const size_t* n, double* ap, size_t* info);
extern void F77_NAME(dpptri)(const char* uplo, const size_t* n, double* ap, size_t* info);
extern void F77_NAME(dpptrs)(const char* uplo, const size_t* n,
		 const size_t* nrhs, const double* ap,
    		 double* b, const size_t* ldb, size_t* info);

extern void F77_NAME(dptcon)(const size_t* n,
		 const double* d, const double* e,
    		 const double* anorm, double* rcond,
    		 double* work, size_t* info);
extern void F77_NAME(dpteqr)(const char* compz, const size_t* n, double* d,
		 double* e, double* z, const size_t* ldz,
    		 double* work, size_t* info);
extern void F77_NAME(dptrfs)(const size_t* n, const size_t* nrhs,
    		 const double* d, const double* e,
    		 const double* df, const double* ef,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx,
    		 double* ferr, double* berr,
    		 double* work, size_t* info);
extern void F77_NAME(dptsv)(const size_t* n, const size_t* nrhs, double* d,
    		double* e, double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dptsvx)(const size_t* fact, const size_t* n,
    		 const size_t* nrhs,
    		 const double* d, const double* e,
    		 double* df, double* ef,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx, double* rcond,
    		 double* ferr, double* berr,
    		 double* work, size_t* info);
extern void F77_NAME(dpttrf)(const size_t* n, double* d, double* e, size_t* info);
extern void F77_NAME(dpttrs)(const size_t* n, const size_t* nrhs,
    		 const double* d, const double* e,
    		 double* b, const size_t* ldb, size_t* info); 
extern void F77_NAME(drscl)(const size_t* n, const double* da,
    		double* x, const size_t* incx);


extern void F77_NAME(dsbev)(const char* jobz, const char* uplo,
    		const size_t* n, const size_t* kd,
    		double* ab, const size_t* ldab,
    		double* w, double* z, const size_t* ldz,
    		double* work, size_t* info);
extern void F77_NAME(dsbevd)(const char* jobz, const char* uplo,
    		 const size_t* n, const size_t* kd,
    		 double* ab, const size_t* ldab,
    		 double* w, double* z, const size_t* ldz,
    		 double* work, const size_t* lwork,
    		 size_t* iwork, const size_t* liwork, size_t* info);
extern void F77_NAME(dsbevx)(const char* jobz, const char* range,
    		 const char* uplo, const size_t* n, const size_t* kd,
    		 double* ab, const size_t* ldab,
    		 double* q, const size_t* ldq,
    		 const double* vl, const double* vu,
    		 const size_t* il, const size_t* iu,
    		 const double* abstol,
    		 size_t* m, double* w,
    		 double* z, const size_t* ldz,
    		 double* work, size_t* iwork,
    		 size_t* ifail, size_t* info);
extern void F77_NAME(dsbgst)(const char* vect, const char* uplo,
    		 const size_t* n, const size_t* ka, const size_t* kb,
    		 double* ab, const size_t* ldab,
    		 double* bb, const size_t* ldbb,
    		 double* x, const size_t* ldx,
    		 double* work, size_t* info);
extern void F77_NAME(dsbgv)(const char* jobz, const char* uplo,
    		const size_t* n, const size_t* ka, const size_t* kb,
    		double* ab, const size_t* ldab,
    		double* bb, const size_t* ldbb,
    		double* w, double* z, const size_t* ldz,
    		double* work, size_t* info);
extern void F77_NAME(dsbtrd)(const char* vect, const char* uplo,
    		 const size_t* n, const size_t* kd,
    		 double* ab, const size_t* ldab,
    		 double* d, double* e,
    		 double* q, const size_t* ldq,
    		 double* work, size_t* info);

extern void F77_NAME(dspcon)(const char* uplo, const size_t* n,
    		 const double* ap, const size_t* ipiv,
    		 const double* anorm, double* rcond,
    		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dspev)(const char* jobz, const char* uplo, const size_t* n,
    		double* ap, double* w, double* z, const size_t* ldz,
    		double* work, size_t* info);
extern void F77_NAME(dspevd)(const char* jobz, const char* uplo,
    		 const size_t* n, double* ap, double* w,
    		 double* z, const size_t* ldz,
    		 double* work, const size_t* lwork,
    		 size_t* iwork, const size_t* liwork, size_t* info);
extern void F77_NAME(dspevx)(const char* jobz, const char* range,
    		 const char* uplo, const size_t* n, double* ap,
    		 const double* vl, const double* vu,
    		 const size_t* il, const size_t* iu,
    		 const double* abstol,
    		 size_t* m, double* w,
    		 double* z, const size_t* ldz,
    		 double* work, size_t* iwork,
    		 size_t* ifail, size_t* info);
extern void F77_NAME(dspgst)(const size_t* itype, const char* uplo,
    		 const size_t* n, double* ap, double* bp, size_t* info);
extern void F77_NAME(dspgv)(const size_t* itype, const char* jobz,
    		const char* uplo, const size_t* n,
    		double* ap, double* bp, double* w,
    		double* z, const size_t* ldz,
    		double* work, size_t* info);

extern void F77_NAME(dsprfs)(const char* uplo, const size_t* n,
    		 const size_t* nrhs, const double* ap,
    		 const double* afp, const size_t* ipiv,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx,
    		 double* ferr, double* berr,
    		 double* work, size_t* iwork, size_t* info);

extern void F77_NAME(dspsv)(const char* uplo, const size_t* n,
    		const size_t* nrhs, double* ap, size_t* ipiv,
    		double* b, const size_t* ldb, size_t* info);

extern void F77_NAME(dspsvx)(const size_t* fact, const char* uplo,
    		 const size_t* n, const size_t* nrhs,
    		 const double* ap, double* afp, size_t* ipiv,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx,
    		 double* rcond, double* ferr, double* berr,
    		 double* work, size_t* iwork, size_t* info);

extern void F77_NAME(dsptrd)(const char* uplo, const size_t* n,
    		 double* ap, double* d, double* e,
    		 double* tau, size_t* info);

extern void F77_NAME(dsptrf)(const char* uplo, const size_t* n,
    		 double* ap, size_t* ipiv, size_t* info);

extern void F77_NAME(dsptri)(const char* uplo, const size_t* n,
    		 double* ap, const size_t* ipiv,
    		 double* work, size_t* info);

extern void F77_NAME(dsptrs)(const char* uplo, const size_t* n,
    		 const size_t* nrhs, const double* ap,
    		 const size_t* ipiv, double* b, const size_t* ldb, size_t* info);


extern void F77_NAME(dstebz)(const char* range, const char* order, const size_t* n,
    		 const double* vl, const double* vu,
    		 const size_t* il, const size_t* iu,
    		 const double *abstol,
    		 const double* d, const double* e,
    		 size_t* m, size_t* nsplit, double* w,
    		 size_t* iblock, size_t* isplit,
    		 double* work, size_t* iwork,
    		 size_t* info);
extern void F77_NAME(dstedc)(const char* compz, const size_t* n,
    		 double* d, double* e,
    		 double* z, const size_t* ldz,
    		 double* work, const size_t* lwork,
    		 size_t* iwork, const size_t* liwork, size_t* info);
extern void F77_NAME(dstein)(const size_t* n, const double* d, const double* e,
    		 const size_t* m, const double* w,
    		 const size_t* iblock, const size_t* isplit,
    		 double* z, const size_t* ldz,
    		 double* work, size_t* iwork,
    		 size_t* ifail, size_t* info);
extern void F77_NAME(dsteqr)(const char* compz, const size_t* n, double* d, double* e,
		 double* z, const size_t* ldz, double* work, size_t* info);
extern void F77_NAME(dsterf)(const size_t* n, double* d, double* e, size_t* info);
extern void F77_NAME(dstev)(const char* jobz, const size_t* n,
    		double* d, double* e,
    		double* z, const size_t* ldz,
    		double* work, size_t* info);
extern void F77_NAME(dstevd)(const char* jobz, const size_t* n,
    		 double* d, double* e,
    		 double* z, const size_t* ldz,
    		 double* work, const size_t* lwork,
    		 size_t* iwork, const size_t* liwork, size_t* info);
extern void F77_NAME(dstevx)(const char* jobz, const char* range,
    		 const size_t* n, double* d, double* e,
    		 const double* vl, const double* vu,
    		 const size_t* il, const size_t* iu,
    		 const double* abstol,
    		 size_t* m, double* w,
    		 double* z, const size_t* ldz,
    		 double* work, size_t* iwork,
    		 size_t* ifail, size_t* info);

extern void F77_NAME(dsycon)(const char* uplo, const size_t* n,
    		 const double* a, const size_t* lda,
    		 const size_t* ipiv,
    		 const double* anorm, double* rcond,
    		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dsyev)(const char* jobz, const char* uplo,
    		const size_t* n, double* a, const size_t* lda,
    		double* w, double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dsyevd)(const char* jobz, const char* uplo,
    		 const size_t* n, double* a, const size_t* lda,
    		 double* w, double* work, const size_t* lwork,
    		 size_t* iwork, const size_t* liwork, size_t* info);
extern void F77_NAME(dsyevx)(const char* jobz, const char* range,
    		 const char* uplo, const size_t* n,
    		 double* a, const size_t* lda,
    		 const double* vl, const double* vu,
    		 const size_t* il, const size_t* iu,
    		 const double* abstol,
    		 size_t* m, double* w,
    		 double* z, const size_t* ldz,
    		 double* work, const size_t* lwork, size_t* iwork,
		 size_t* ifail, size_t* info);
extern void F77_NAME(dsyevr)(const char *jobz, const char *range, const char *uplo,
		 const size_t *n, double *a, const size_t *lda,
		 const double *vl, const double *vu,
		 const size_t *il, const size_t *iu,
		 const double *abstol, size_t *m, double *w, 
		 double *z, const size_t *ldz, size_t *isuppz, 
		 double *work, const size_t *lwork,
		 size_t *iwork, const size_t *liwork,
		 size_t *info);
extern void F77_NAME(dsygs2)(const size_t* itype, const char* uplo,
    		 const size_t* n, double* a, const size_t* lda,
    		 const double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dsygst)(const size_t* itype, const char* uplo,
    		 const size_t* n, double* a, const size_t* lda,
    		 const double* b, const size_t* ldb, size_t* info);
extern void F77_NAME(dsygv)(const size_t* itype, const char* jobz,
    		const char* uplo, const size_t* n,
    		double* a, const size_t* lda,
    		double* b, const size_t* ldb,
    		double* w, double* work, const size_t* lwork,
    		size_t* info);
extern void F77_NAME(dsyrfs)(const char* uplo, const size_t* n,
    		 const size_t* nrhs,
    		 const double* a, const size_t* lda,
    		 const double* af, const size_t* ldaf,
    		 const size_t* ipiv,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx,
    		 double* ferr, double* berr,
    		 double* work, size_t* iwork, size_t* info);

extern void F77_NAME(dsysv)(const char* uplo, const size_t* n,
    		const size_t* nrhs,
    		double* a, const size_t* lda, size_t* ipiv,
    		double* b, const size_t* ldb,
    		double* work, const size_t* lwork, size_t* info);

extern void F77_NAME(dsysvx)(const size_t* fact, const char* uplo,
    		 const size_t* n, const size_t* nrhs,
    		 const double* a, const size_t* lda,
    		 double* af, const size_t* ldaf, size_t* ipiv,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx, double* rcond,
    		 double* ferr, double* berr,
    		 double* work, const size_t* lwork,
    		 size_t* iwork, size_t* info);

extern void F77_NAME(dsytd2)(const char* uplo, const size_t* n,
    		 double* a, const size_t* lda,
    		 double* d, double* e, double* tau,
    		 size_t* info);

extern void F77_NAME(dsytf2)(const char* uplo, const size_t* n,
    		 double* a, const size_t* lda,
    		 size_t* ipiv, size_t* info);

extern void F77_NAME(dsytrd)(const char* uplo, const size_t* n,
    		 double* a, const size_t* lda,
    		 double* d, double* e, double* tau,
    		 double* work, const size_t* lwork, size_t* info);

extern void F77_NAME(dsytrf)(const char* uplo, const size_t* n,
    		 double* a, const size_t* lda, size_t* ipiv,
    		 double* work, const size_t* lwork, size_t* info);

extern void F77_NAME(dsytri)(const char* uplo, const size_t* n,
    		 double* a, const size_t* lda, const size_t* ipiv,
    		 double* work, size_t* info); 

extern void F77_NAME(dsytrs)(const char* uplo, const size_t* n,
    		 const size_t* nrhs,
    		 const double* a, const size_t* lda,
    		 const size_t* ipiv,
    		 double* b, const size_t* ldb, size_t* info);

extern void F77_NAME(dtbcon)(const char* norm, const char* uplo,
    		 const char* diag, const size_t* n, const size_t* kd,
    		 const double* ab, const size_t* ldab,
    		 double* rcond, double* work,
    		 size_t* iwork, size_t* info); 
extern void F77_NAME(dtbrfs)(const char* uplo, const char* trans,
    		 const char* diag, const size_t* n, const size_t* kd,
    		 const size_t* nrhs,
    		 const double* ab, const size_t* ldab,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx,
    		 double* ferr, double* berr,
    		 double* work, size_t* iwork, size_t* info);  
extern void F77_NAME(dtbtrs)(const char* uplo, const char* trans,
    		 const char* diag, const size_t* n,
    		 const size_t* kd, const size_t* nrhs,
    		 const double* ab, const size_t* ldab,
    		 double* b, const size_t* ldb, size_t* info); 

extern void F77_NAME(dtgevc)(const char* side, const char* howmny,
    		 const size_t* select, const size_t* n,
    		 const double* a, const size_t* lda,
    		 const double* b, const size_t* ldb,
    		 double* vl, const size_t* ldvl,
    		 double* vr, const size_t* ldvr,
    		 const size_t* mm, size_t* m, double* work, size_t* info);

extern void F77_NAME(dtgsja)(const char* jobu, const char* jobv, const char* jobq,
    		 const size_t* m, const size_t* p, const size_t* n,
    		 const size_t* k, const size_t* l,
    		 double* a, const size_t* lda,
    		 double* b, const size_t* ldb,
    		 const double* tola, const double* tolb,
    		 double* alpha, double* beta,
    		 double* u, const size_t* ldu,
    		 double* v, const size_t* ldv,
    		 double* q, const size_t* ldq,
    		 double* work, size_t* ncycle, size_t* info);
extern void F77_NAME(dtpcon)(const char* norm, const char* uplo,
    		 const char* diag, const size_t* n,
    		 const double* ap, double* rcond,
    		 double* work, size_t* iwork, size_t* info);

extern void F77_NAME(dtprfs)(const char* uplo, const char* trans,
    		 const char* diag, const size_t* n,
    		 const size_t* nrhs, const double* ap,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx,
    		 double* ferr, double* berr,
    		 double* work, size_t* iwork, size_t* info);

extern void F77_NAME(dtptri)(const char* uplo, const char* diag,
    		 const size_t* n, double* ap, size_t* info);

extern void F77_NAME(dtptrs)(const char* uplo, const char* trans,
    		 const char* diag, const size_t* n,
    		 const size_t* nrhs, const double* ap,
    		 double* b, const size_t* ldb, size_t* info); 

extern void F77_NAME(dtrcon)(const char* norm, const char* uplo,
    		 const char* diag, const size_t* n,
    		 const double* a, const size_t* lda,
    		 double* rcond, double* work,
    		 size_t* iwork, size_t* info);

extern void F77_NAME(dtrevc)(const char* side, const char* howmny,
    		 const size_t* select, const size_t* n,
    		 const double* t, const size_t* ldt,
    		 double* vl, const size_t* ldvl,
    		 double* vr, const size_t* ldvr,
    		 const size_t* mm, size_t* m, double* work, size_t* info);

extern void F77_NAME(dtrexc)(const char* compq, const size_t* n,
    		 double* t, const size_t* ldt,
    		 double* q, const size_t* ldq,
    		 size_t* ifst, size_t* ILST,
    		 double* work, size_t* info);

extern void F77_NAME(dtrrfs)(const char* uplo, const char* trans,
    		 const char* diag, const size_t* n, const size_t* nrhs,
    		 const double* a, const size_t* lda,
    		 const double* b, const size_t* ldb,
    		 double* x, const size_t* ldx,
    		 double* ferr, double* berr,
    		 double* work, size_t* iwork, size_t* info); 

extern void F77_NAME(dtrsen)(const char* job, const char* compq,
    		 const size_t* select, const size_t* n,
    		 double* t, const size_t* ldt,
    		 double* q, const size_t* ldq,
    		 double* wr, double* wi,
    		 size_t* m, double* s, double* sep,
    		 double* work, const size_t* lwork,
    		 size_t* iwork, const size_t* liwork, size_t* info);

extern void F77_NAME(dtrsna)(const char* job, const char* howmny,
    		 const size_t* select, const size_t* n,
    		 const double* t, const size_t* ldt,
    		 const double* vl, const size_t* ldvl,
    		 const double* vr, const size_t* ldvr,
    		 double* s, double* sep, const size_t* mm,
    		 size_t* m, double* work, const size_t* lwork,
    		 size_t* iwork, size_t* info);

extern void F77_NAME(dtrsyl)(const char* trana, const char* tranb,
    		 const size_t* isgn, const size_t* m, const size_t* n,
    		 const double* a, const size_t* lda,
    		 const double* b, const size_t* ldb,
    		 double* c, const size_t* ldc,
    		 double* scale, size_t* info);  

extern void F77_NAME(dtrti2)(const char* uplo, const char* diag,
    		 const size_t* n, double* a, const size_t* lda,
    		 size_t* info); 

extern void F77_NAME(dtrtri)(const char* uplo, const char* diag,
    		 const size_t* n, double* a, const size_t* lda,
    		 size_t* info); 

extern void F77_NAME(dtrtrs)(const char* uplo, const char* trans,
    		 const char* diag, const size_t* n, const size_t* nrhs,
    		 const double* a, const size_t* lda,
    		 double* b, const size_t* ldb, size_t* info); 

extern void F77_NAME(dtzrqf)(const size_t* m, const size_t* n,
    		 double* a, const size_t* lda,
    		 double* tau, size_t* info);



extern void F77_NAME(dhgeqz)(const char* job, const char* compq, const char* compz,
		 const size_t* n, const size_t *ILO, const size_t* IHI,
		 double* a, const size_t* lda,
		 double* b, const size_t* ldb,
		 double* alphar, double* alphai, const double* beta,
		 double* q, const size_t* ldq,
		 double* z, const size_t* ldz,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dhsein)(const char* side, const char* eigsrc,
		 const char* initv, size_t* select,
		 const size_t* n, double* h, const size_t* ldh,
		 double* wr, double* wi,
		 double* vl, const size_t* ldvl,
		 double* vr, const size_t* ldvr,
		 const size_t* mm, size_t* m, double* work,
		 size_t* ifaill, size_t* ifailr, size_t* info);
extern void F77_NAME(dhseqr)(const char* job, const char* compz, const size_t* n,
		 const size_t* ilo, const size_t* ihi,
		 double* h, const size_t* ldh,
		 double* wr, double* wi,
		 double* z, const size_t* ldz,
		 double* work, const size_t* lwork, size_t* info);
extern void F77_NAME(dlabad)(double* small, double* large);
extern void F77_NAME(dlabrd)(const size_t* m, const size_t* n, const size_t* nb,
		 double* a, const size_t* lda, double* d, double* e,
		 double* tauq, double* taup,
		 double* x, const size_t* ldx, double* y, const size_t* ldy);
extern void F77_NAME(dlacon)(const size_t* n, double* v, double* x,
		 size_t* isgn, double* est, size_t* kase);
extern void F77_NAME(dlacpy)(const char* uplo, const size_t* m, const size_t* n,
		 const double* a, const size_t* lda,
		 double* b, const size_t* ldb);
extern void F77_NAME(dladiv)(const double* a, const double* b,
		 const double* c, const double* d,
		 double* p, double* q);
extern void F77_NAME(dlae2)(const double* a, const double* b, const double* c,
		double* rt1, double* rt2);
extern void F77_NAME(dlaebz)(const size_t* ijob, const size_t* nitmax, const size_t* n,
		 const size_t* mmax, const size_t* minp, const size_t* nbmin,
		 const double* abstol, const double* reltol,
		 const double* pivmin, double* d, double* e,
		 double* e2, size_t* nval, double* ab, double* c,
		 size_t* mout, size_t* nab, double* work, size_t* iwork,
		 size_t* info);
extern void F77_NAME(dlaed0)(const size_t* icompq, const size_t* qsiz, const size_t* n,
		 double* d, double* e, double* q, const size_t* ldq,
		 double* qstore, const size_t* ldqs,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dlaed1)(const size_t* n, double* d, double* q, const size_t* ldq,
		 size_t* indxq, const double* rho, const size_t* cutpnt,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dlaed2)(const size_t* k, const size_t* n, double* d,
		 double* q, const size_t* ldq, size_t* indxq,
		 double* rho, const size_t* cutpnt, double* z,
		 double* dlamda, double* q2, const size_t *ldq2,
		 size_t* indxc, size_t* w, size_t* indxp, size_t* indx,
		 size_t* coltyp, size_t* info);
extern void F77_NAME(dlaed3)(const size_t* k, const size_t* kstart,
		 const size_t *kstop, const size_t* n,
		 double* d, double* q, const size_t* ldq,
		 const double* rho, const size_t* cutpnt,
		 double* dlamda, size_t* q2, const size_t* ldq2,
		 size_t* indxc, size_t* ctot, double* w,
		 double* s, const size_t* lds, size_t* info);
extern void F77_NAME(dlaed4)(const size_t* n, const size_t* i, const double* d,
		 const double* z, const double* delta,
		 const double* rho, double* dlam, size_t* info);
extern void F77_NAME(dlaed5)(const size_t* i, const double* d, const double* z,
		 double* delta, const double* rho, double* dlam);
extern void F77_NAME(dlaed6)(const size_t* kniter, const size_t* orgati,
		 const double* rho, const double* d,
		 const double* z, const double* finit,
		 double* tau, size_t* info);
extern void F77_NAME(dlaed7)(const size_t* icompq, const size_t* n,
		 const size_t* qsiz, const size_t* tlvls,
		 const size_t* curlvl, const size_t* curpbm,
		 double* d, double* q, const size_t* ldq,
		 size_t* indxq, const double* rho, const size_t* cutpnt,
		 double* qstore, double* qptr, const size_t* prmptr,
		 const size_t* perm, const size_t* givptr,
		 const size_t* givcol, const double* givnum,
		 double* work, size_t* iwork, size_t* info);
extern void F77_NAME(dlaed8)(const size_t* icompq, const size_t* k,
		 const size_t* n, const size_t* qsiz,
		 double* d, double* q, const size_t* ldq,
		 const size_t* indxq, double* rho,
		 const size_t* cutpnt, const double* z,
		 double* dlamda, double* q2, const size_t* ldq2,
		 double* w, size_t* perm, size_t* givptr,
		 size_t* givcol, double* givnum, size_t* indxp,
		 size_t* indx, size_t* info);
extern void F77_NAME(dlaed9)(const size_t* k, const size_t* kstart, const size_t* kstop,
		 const size_t* n, double* d, double* q, const size_t* ldq,
		 const double* rho, const double* dlamda,
		 const double* w, double* s, const size_t* lds, size_t* info);
extern void F77_NAME(dlaeda)(const size_t* n, const size_t* tlvls, const size_t* curlvl,
		 const size_t* curpbm, const size_t* prmptr, const size_t* perm,
		 const size_t* givptr, const size_t* givcol,
		 const double* givnum, const double* q,
		 const size_t* qptr, double* z, double* ztemp, size_t* info);
extern void F77_NAME(dlaein)(const size_t* rightv, const size_t* noinit, const size_t* n,
		 const double* h, const size_t* ldh,
		 const double* wr, const double* wi,
		 double* vr, double* vi,
		 double* b, const size_t* ldb, double* work,
		 const double* eps3, const double* smlnum,
		 const double* bignum, size_t* info);
extern void F77_NAME(dlaev2)(const double* a, const double* b, const double* c, 
		 double* rt1, double* rt2, double* cs1, double *sn1);
extern void F77_NAME(dlaexc)(const size_t* wantq, const size_t* n, double* t, const size_t* ldt,
    		  double* q, const size_t* ldq, const size_t* j1,
		 const size_t* n1, const size_t* n2, double* work, size_t* info);
extern void F77_NAME(dlag2)(const double* a, const size_t* lda, const double* b,
		const size_t* ldb, const double* safmin,
		double* scale1, double* scale2,
		double* wr1, double* wr2, double* wi);
extern void F77_NAME(dlags2)(const size_t* upper,
		 const double* a1, const double* a2, const double* a3,
		 const double* b1, const double* b2, const double* b3,
		 double* csu, double* snu,
		 double* csv, double* snv, double *csq, double *snq);
extern void F77_NAME(dlagtf)(const size_t* n, double* a, const double* lambda,
		 double* b, double* c, const double *tol,
		 double* d, size_t* in, size_t* info);
extern void F77_NAME(dlagtm)(const char* trans, const size_t* n, const size_t* nrhs,
		 const double* alpha, const double* dl,
		 const double* d, const double* du,
		 const double* x, const size_t* ldx, const double* beta,
		 double* b, const size_t* ldb);
extern void F77_NAME(dlagts)(const size_t* job, const size_t* n,
		 const double* a, const double* b,
		 const double* c, const double* d,
		 const size_t* in, double* y, double* tol, size_t* info);
extern void F77_NAME(dlahqr)(const size_t* wantt, const size_t* wantz, const size_t* n,
		 const size_t* ilo, const size_t* ihi,
		 double* H, const size_t* ldh, double* wr, double* wi,
		 const size_t* iloz, const size_t* ihiz,
		 double* z, const size_t* ldz, size_t* info);
extern void F77_NAME(dlahrd)(const size_t* n, const size_t* k, const size_t* nb,
		 double* a, const size_t* lda,
		 double* tau, double* t, const size_t* ldt,
		 double* y, const size_t* ldy);
extern void F77_NAME(dlaic1)(const size_t* job, const size_t* j, const double* x,
		 const double* sest, const double* w,
		 const double* gamma, double* sestpr,
		 double* s, double* c);
extern void F77_NAME(dlaln2)(const size_t* ltrans, const size_t* na, const size_t* nw,
		 const double* smin, const double* ca,
		 const double* a, const size_t* lda,
		 const double* d1, const double* d2,
		 const double* b, const size_t* ldb,
		 const double* wr, const double* wi,
		 double* x, const size_t* ldx, double* scale,
		 double* xnorm, size_t* info);
extern double F77_NAME(dlamch)(const char* cmach);
extern void F77_NAME(dlamrg)(const size_t* n1, const size_t* n2, const double* a,
		 const size_t* dtrd1, const size_t* dtrd2, size_t* index);
extern double F77_NAME(dlangb)(const char* norm, const size_t* n,
		 const size_t* kl, const size_t* ku, const double* ab,
		 const size_t* ldab, double* work);
extern double F77_NAME(dlange)(const char* norm, const size_t* m, const size_t* n,
		 const double* a, const size_t* lda, double* work);
extern double F77_NAME(dlangt)(const char* norm, const size_t* n,
		 const double* dl, const double* d,
		 const double* du);
extern double F77_NAME(dlanhs)(const char* norm, const size_t* n,
		 const double* a, const size_t* lda, double* work);
extern double F77_NAME(dlansb)(const char* norm, const char* uplo,
		 const size_t* n, const size_t* k,
		 const double* ab, const size_t* ldab, double* work);
extern double F77_NAME(dlansp)(const char* norm, const char* uplo,
		 const size_t* n, const double* ap, double* work);
extern double F77_NAME(dlanst)(const char* norm, const size_t* n,
		 const double* d, const double* e);
extern double F77_NAME(dlansy)(const char* norm, const char* uplo, const size_t* n,
		 const double* a, const size_t* lda, double* work);
extern double F77_NAME(dlantb)(const char* norm, const char* uplo,
		 const char* diag, const size_t* n, const size_t* k,
		 const double* ab, const size_t* ldab, double* work);
extern double F77_NAME(dlantp)(const char* norm, const char* uplo, const char* diag,
		 const size_t* n, const double* ap, double* work);
extern double F77_NAME(dlantr)(const char* norm, const char* uplo,
		 const char* diag, const size_t* m, const size_t* n,
		 const double* a, const size_t* lda, double* work);
extern void F77_NAME(dlanv2)(double* a, double* b, double* c, double* d,
		 double* rt1r, double* rt1i, double* rt2r, double* rt2i,
		 double* cs, double *sn);
extern void F77_NAME(dlapll)(const size_t* n, double* x, const size_t* incx,
		 double* y, const size_t* incy, double* ssmin);
extern void F77_NAME(dlapmt)(const size_t* forwrd, const size_t* m, const size_t* n,
		 double* x, const size_t* ldx, const size_t* k);
extern double F77_NAME(dlapy2)(const double* x, const double* y);
extern double F77_NAME(dlapy3)(const double* x, const double* y, const double* z);
extern void F77_NAME(dlaqgb)(const size_t* m, const size_t* n,
		 const size_t* kl, const size_t* ku,
		 double* ab, const size_t* ldab,
		 double* r, double* c,
		 double* rowcnd, double* colcnd,
		 const double* amax, char* equed);
extern void F77_NAME(dlaqge)(const size_t* m, const size_t* n,
		 double* a, const size_t* lda,
		 double* r, double* c,
		 double* rowcnd, double* colcnd,
		 const double* amax, char* equed);
extern void F77_NAME(dlaqsb)(const char* uplo, const size_t* n, const size_t* kd,
		 double* ab, const size_t* ldab, const double* s,
		 const double* scond, const double* amax, char* equed);
extern void F77_NAME(dlaqsp)(const char* uplo, const size_t* n,
		 double* ap, const double* s, const double* scond,
		 const double* amax, size_t* equed);
extern void F77_NAME(dlaqsy)(const char* uplo, const size_t* n,
		 double* a, const size_t* lda,
		 const double* s, const double* scond, 
		 const double* amax, size_t* equed);
extern void F77_NAME(dlaqtr)(const size_t* ltran, const size_t* lreal, const size_t* n,
		 const double* t, const size_t* ldt,
		 const double* b, const double* w,
		 double* scale, double* x, double* work, size_t* info);
extern void F77_NAME(dlar2v)(const size_t* n, double* x, double* y,
		 double* z, const size_t* incx,
		 const double* c, const double* s,
		 const size_t* incc);
extern void F77_NAME(dlarf)(const char* side, const size_t* m, const size_t* n,
		const double* v, const size_t* incv, const double* tau,
		double* c, const size_t* ldc, double* work);
extern void F77_NAME(dlarfb)(const char* side, const char* trans,
		 const char* direct, const char* storev,
		 const size_t* m, const size_t* n, const size_t* k,
		 const double* v, const size_t* ldv,
		 const double* t, const size_t* ldt,
		 double* c, const size_t* ldc,
		 double* work, const size_t* lwork);
extern void F77_NAME(dlarfg)(const size_t* n, const double* alpha,
		 double* x, const size_t* incx, double* tau);
extern void F77_NAME(dlarft)(const char* direct, const char* storev,
		 const size_t* n, const size_t* k, double* v, const size_t* ldv,
		 const double* tau, double* t, const size_t* ldt);
extern void F77_NAME(dlarfx)(const char* side, const size_t* m, const size_t* n,
		 const double* v, const double* tau,
		 double* c, const size_t* ldc, double* work);
extern void F77_NAME(dlargv)(const size_t* n, double* x, const size_t* incx,
		 double* y, const size_t* incy, double* c, const size_t* incc);
extern void F77_NAME(dlarnv)(const size_t* idist, size_t* iseed, const size_t* n, double* x);
extern void F77_NAME(dlartg)(const double* f, const double* g, double* cs,
		 double* sn, double *r);
extern void F77_NAME(dlartv)(const size_t* n, double* x, const size_t* incx,
		 double* y, const size_t* incy,
		 const double* c, const double* s,
		 const size_t* incc);
extern void F77_NAME(dlaruv)(size_t* iseed, const size_t* n, double* x);

extern void F77_NAME(dlas2)(const double* f, const double* g, const double* h,
    		 double* ssmin, double* ssmax);

extern void F77_NAME(dlascl)(const char* type,
		 const size_t* kl,const size_t* ku,
		 double* cfrom, double* cto,
		 const size_t* m, const size_t* n,
		 double* a, const size_t* lda, size_t* info);
    
extern void F77_NAME(dlaset)(const char* uplo, const size_t* m, const size_t* n,
		 const double* alpha, const double* beta,
		 double* a, const size_t* lda);
extern void F77_NAME(dlasq1)(const size_t* n, double* d, double* e,
		 double* work, size_t* info);
extern void F77_NAME(dlasq2)(const size_t* m, double* q, double* e,
		 double* qq, double* ee, const double* eps,
		 const double* tol2, const double* small2,
		 double* sup, size_t* kend, size_t* info);
extern void F77_NAME(dlasq3)(size_t* n, double* q, double* e, double* qq,
		 double* ee, double* sup, double *sigma,
		 size_t* kend, size_t* off, size_t* iphase,
		 const size_t* iconv, const double* eps,
		 const double* tol2, const double* small2);
extern void F77_NAME(dlasq4)(const size_t* n, const double* q, const double* e,
		 double* tau, double* sup);
extern void F77_NAME(dlasr)(const char* side, const char* pivot,
		const char* direct, const size_t* m, const size_t* n,
		const double* c, const double* s,
		double* a, const size_t* lda);
extern void F77_NAME(dlasrt)(const char* id, const size_t* n, double* d, size_t* info);
extern void F77_NAME(dlassq)(const size_t* n, const double* x, const size_t* incx,
		 double* scale, double* sumsq);
extern void F77_NAME(dlasv2)(const double* f, const double* g, const double* h,
		 double* ssmin, double* ssmax, double* snr, double* csr,
		 double* snl, double* csl);
extern void F77_NAME(dlaswp)(const size_t* n, double* a, const size_t* lda,
		 const size_t* k1, const size_t* k2,
		 const size_t* ipiv, const size_t* incx);
extern void F77_NAME(dlasy2)(const size_t* ltranl, const size_t* ltranr,
		 const size_t* isgn, const size_t* n1, const size_t* n2,
		 const double* tl, const size_t* ldtl,
		 const double* tr, const size_t* ldtr,
		 const double* b, const size_t* ldb,
		 double* scale, double* x, const size_t* ldx,
		 double* xnorm, size_t* info);
extern void F77_NAME(dlasyf)(const char* uplo, const size_t* n,
		 const size_t* nb, const size_t* kb,
		 double* a, const size_t* lda, size_t* ipiv,
		 double* w, const size_t* ldw, size_t* info);
extern void F77_NAME(dlatbs)(const char* uplo, const char* trans,
		 const char* diag, const char* normin,
		 const size_t* n, const size_t* kd,
		 const double* ab, const size_t* ldab,
		 double* x, double* scale, double* cnorm, size_t* info);
extern void F77_NAME(dlatps)(const char* uplo, const char* trans,
		 const char* diag, const char* normin,
		 const size_t* n, const double* ap,
		 double* x, double* scale, double* cnorm, size_t* info);
extern void F77_NAME(dlatrd)(const char* uplo, const size_t* n, const size_t* nb,
		 double* a, const size_t* lda, double* e, double* tau,
		 double* w, const size_t* ldw);
extern void F77_NAME(dlatrs)(const char* uplo, const char* trans,
		 const char* diag, const char* normin,
		 const size_t* n, const double* a, const size_t* lda,
		 double* x, double* scale, double* cnorm, size_t* info);
extern void F77_NAME(dlatzm)(const char* side, const size_t* m, const size_t* n,
		 const double* v, const size_t* incv,
		 const double* tau, double* c1, double* c2,
		 const size_t* ldc, double* work);
extern void F77_NAME(dlauu2)(const char* uplo, const size_t* n,
		 double* a, const size_t* lda, size_t* info);
extern void F77_NAME(dlauum)(const char* uplo, const size_t* n,
		 double* a, const size_t* lda, size_t* info);


#ifdef  __cplusplus
}
#endif

#endif /* LAPACK_H */
