#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void train_regression_stump_R(const int *N_train, const int *P_test,
		const double *x_train, const double *y_train, const double *s_f,
		const double *s_s, const int *lf, const int *M_train, const double *nu, double *I,
		double *f0, int *featI, double * featT, double *gamma_l,
		double *gamma_r);

void test_regression_stump_R(const int *N_test, const int *P_test,
		const int *P_train, const double *x_test, const double *y_test,
		const int *M_test, const int *M_train, const double *nu,
		const double *f0, const int *featI, const double * featT,
		const double *gamma_l, const double *gamma_r, double *loss, double *p);

static R_NativePrimitiveArgType train_regression_stump_R_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, 
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType test_regression_stump_R_t[] = {
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_CMethodDef cMethods[] = {

  {"train_regression_stump_R", (DL_FUNC) &train_regression_stump_R, 15, train_regression_stump_R_t},
  {"test_regression_stump_R", (DL_FUNC) &test_regression_stump_R, 15, test_regression_stump_R_t}, 
  {NULL,NULL,0}
  
};

void R_init_RGBM(DllInfo *info)
{
  /* Register routines,
     allocate resources. */
	R_registerRoutines(info,cMethods,NULL,NULL,NULL);
  R_useDynamicSymbols(info,FALSE);
}

