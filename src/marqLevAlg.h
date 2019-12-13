#include <R_ext/RS.h>

void F77_SUB(dsinv)(double * a,
		    int * n,
		    double * eps,
		    int * ier,
		    double * det);

void F77_SUB(dchole)(double * a,
		     int * k,
		     int * nq,
		     int * idpos);
