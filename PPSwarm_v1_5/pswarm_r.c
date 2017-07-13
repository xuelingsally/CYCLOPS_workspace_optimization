#include "pswarm_r.h"

static jmp_buf Jb;


void catchfpe(int n)
{
  Rprintf("\nFloating point error.\n");
  fflush(stdout);
  longjmp(Jb,1);
}


void dim_check(char *txt, SEXP obj, int x, int y)
{
  SEXP r_dim;
  int nx, ny;

  if(obj==R_NilValue)
    return;

  r_dim = getAttrib(obj, R_DimSymbol);
  if(r_dim==R_NilValue)
    return;
    
  nx=INTEGER(r_dim)[0];
  ny=INTEGER(r_dim)[1];
    
  if(nx!=x)
    error("%s must have %d row(s) and as %d\n",txt,x, nx);
    
  if(ny!=y)
    error("%s must have %d column(s) and as %d\n", txt, y, ny);

}

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}


SEXP createRRealVector(int n, double* x) {

  SEXP ans;
  int i;
  
  PROTECT(ans = allocVector(REALSXP, n));
  for (i = 0; i < n; i++)
    if(x){
      REAL(ans)[i] = x[i];
    } else {
      REAL(ans)[i] = 0.0;
    }
  UNPROTECT(1);
  
  return ans;
}

SEXP createRRealMatrix(int n, int m, double* x) {

  SEXP ans, dim;
  int i;
  
  PROTECT(ans = allocVector(REALSXP, n*m));

  memcpy(&REAL(ans)[0],x,m*n*sizeof(double));

  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = n;
  INTEGER(dim)[1] = m;
  setAttrib(ans, R_DimSymbol, dim);

  UNPROTECT(2);
  
  return ans;
}


SEXP createRRealScalar(double x) {

  SEXP ans;
  
  PROTECT(ans = allocVector(REALSXP, 1));
  REAL(ans)[0] = x;
  UNPROTECT(1);
  
  return ans;
}


SEXP createRIntScalar(int x) {

  SEXP ans;

  PROTECT(ans = allocVector(INTSXP, 1));
  INTEGER(ans)[0] = x;
  UNPROTECT(1);

  return ans;

}


void r_objfun(int n, int m, double *x, double *lb, double *ub, double *fx)
{
  SEXP r_fun, t, s, res;
  int i,j,k;
  double *xx;

  if(x==NULL || m==0)
    return;

  xx=(double *)malloc(n*m*sizeof(double));
  memcpy(xx,x,n*m*sizeof(double));

  /* pswarm controls bound feasibility, but just in case... */
  for(j=0;j<m;j++){
	for(i=0;i<n;i++){
		if(xx[j*n+i]<lb[i] || xx[j*n+i]>ub[i]){
			Rprintf("Error computing objective function for unfeasible bound point\nReturning all as infinity.");
			for(k=0;k<m;k++)
				fx[k]=1e20;
			free(xx);
			return;
		}
	}
  }

  /* Build a list for calling the objective function */
  PROTECT(t = s = allocList(2));
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, r_objf); t = CDR(t);
  SETCAR(t, createRRealMatrix(n,m,xx));
  PROTECT(r_fun=eval(s, r_environment));
  UNPROTECT(2);

  if(r_fun!=R_UnboundValue){  
	  PROTECT(res = eval(r_fun,r_environment));
	  if(!isVector(res)){
		  Rprintf("Error evaluating objective function!!\nA vector should be returned\n");
		  for(k=0;k<m;k++)
				fx[k]=1e20;
			return;
	  }
	  for(k=0;k<m;k++)
		  fx[k]= REAL(res)[k];
	  UNPROTECT(1);
  } else {
    Rprintf("Error evaluating objective function!!\nLet's hope it's fine!!");
	for(k=0;k<m;k++)
		fx[k]=1e20;
  }

  free(xx);
  return;
}


double r_outfcn(int n, int s, int iter, int gbest, struct swarm *pop)
{
	SEXP r_out, tt, ts;
	double fx;

	if(r_outf && r_outf!=R_NilValue){
	  PROTECT(tt = ts = allocList(5));
	  SET_TYPEOF(ts, LANGSXP);
	  SETCAR(tt, r_outf); tt=CDR(tt);
	  SETCAR(tt, createRIntScalar(iter)); tt=CDR(tt);
	  SETCAR(tt, createRIntScalar(gbest)); tt=CDR(tt);
	  SETCAR(tt, createRRealScalar((pop->fy[gbest]))); tt=CDR(tt);
	  SETCAR(tt, createRRealVector(n,&(pop->y[gbest*n])));
	  PROTECT(r_out=eval(ts, r_environment));
	  UNPROTECT(2);

	  if(r_out!=R_UnboundValue && r_out!=R_NilValue){
	   	fx=(REAL(eval(r_out,r_environment)))[0];
	  } else {
	  	fx= 1.0;
	  }

	 return(fx);

	} else {
		  /* print to stdout */
	  	if(iter==0){
			Rprintf("\n  Iter     Leader     Objective  ");
			Rprintf("\n  -------------------------------\n");
		}

		Rprintf("    %4d   %4d   %4.6e\n", iter, gbest, pop->fy[gbest]);

		return(1.0);
	}

	return(1.0);
}



SEXP pswarm_r(SEXP r_problem_ptr, SEXP r_options_ptr, SEXP environment_ptr)
{
  SEXP r_A, r_b, r_lb, r_ub, r_x0, A_dim, r_n;
  SEXP r_opt, s, t, r_fobj, r_x, r_ret;
  double *lb=NULL, *ub=NULL;
  double *A=NULL, *b=NULL, *x0=NULL;
  int n, i, j, lincons, nx, ny, exit_code;
  double *sol=NULL;
  double f;

  if(!saved_options){
    memcpy(&opt_backup, &opt, sizeof(struct Options));
    saved_options++;
    //    Rprintf("Saving options\n");
  } else {
    memcpy(&opt, &opt_backup, sizeof(struct Options));
    //    Rprintf("Recovering saved options\n");
  }

  /* Some checkup on the input parameters */

  /* Last parameter must be an environment */
  if(!isEnvironment(environment_ptr))
    error("Third parameter should be an environment");

  r_environment=environment_ptr; /* global setting */

  /* First parameter is a list with the problem definition */
  if(!isVector(r_problem_ptr))
    error("First parameter should be a list with the problem definition\n");

  //  PrintValue(r_problem_ptr);

  /* At least an objective function pointer should be provided */
  if((r_objf = getListElement(r_problem_ptr, "objf")) == R_NilValue)
     error("The objective function name must be provided\n");

  /* check for number of variables */
  if((r_n = getListElement(r_problem_ptr, "Variables")) == R_NilValue)
    error("The number of variables must be provided\n");

  n=(INTEGER(AS_INTEGER(r_n)))[0];
  if(n<=0)
    error("Number of variables must be a positive integer\n");

  /* get lower and upper bounds */
  r_lb = getListElement(r_problem_ptr, "lb");
  r_ub = getListElement(r_problem_ptr, "ub");

  dim_check("lb", r_lb, n, 1); /* exits on error */
  dim_check("ub", r_ub, n, 1); /* exits on error */

  if((lb=(double *) malloc(n*sizeof(double)))==NULL)
    error("Out of memory\n");
  if((ub=(double *) malloc(n*sizeof(double)))==NULL)
    error("Out of memory\n");

  if(r_lb == R_NilValue){
    for(i=0;i<n;i++)
      lb[i]=-1e20;
  } else {
    for(i=0;i<n;i++)
      lb[i]=(REAL(r_lb))[i];
  }

  if(r_ub == R_NilValue){
    for(i=0;i<n;i++)
      ub[i]=1e20;
  } else {
    for(i=0;i<n;i++)
      ub[i]=(REAL(r_ub))[i];
  }



  /* get linear constraints */
  r_A = getListElement(r_problem_ptr, "A");
  r_b = getListElement(r_problem_ptr, "b");

  if(r_A == R_NilValue || r_b == R_NilValue){
    lincons=0;
    A=NULL;
    b=NULL;
  } else {
    A_dim = getAttrib(r_A, R_DimSymbol);

    nx=INTEGER(A_dim)[0];
    ny=INTEGER(A_dim)[1];

    if(ny!=n)
      error("List element A must be a matrix with %d columns\n",n);
    
    lincons=nx;

    dim_check("b", r_b, nx, 1); /* exits on error */    

    if((A=(double *) malloc(n*lincons*sizeof(double)))==NULL)
      error("Out of memory\n");
    if((b=(double *) malloc(lincons*sizeof(double)))==NULL)
      error("Out of memory\n");
    
    for(j=0;j<lincons;j++){
      b[j]=REAL(r_b)[j];
      for(i=0;i<n;i++)
	A[i+j*n]=REAL(r_A)[i+j*n];
    }
  }

  /* get initial guess, if provided */
  r_x0 = getListElement(r_problem_ptr, "x0");

  if(r_x0 == R_NilValue){
    x0 = NULL;
  } else {
    dim_check("x0", r_x0, n, 1); /* exits on error */    

    if((x0=(double *) malloc(n*sizeof(double)))==NULL)
      error("Out of memory\n");

    for(i=0;i<n;i++){
      x0[i]=REAL(r_x0)[i];
    }

  }
 
  //PrintValue(r_options_ptr);

  r_opt = getListElement(r_options_ptr, "cognitial");
  if(r_opt != R_NilValue)
    opt.mu=REAL(r_opt)[0];

  r_opt = getListElement(r_options_ptr, "fweight");
  if(r_opt != R_NilValue)
    opt.fweight=REAL(r_opt)[0];

  r_opt = getListElement(r_options_ptr, "iweight");
  if(r_opt != R_NilValue)
    opt.iweight=REAL(r_opt)[0];

  r_opt = getListElement(r_options_ptr, "maxf");
  if(r_opt != R_NilValue)
    opt.maxf=(int)(REAL(r_opt)[0]);

  r_opt = getListElement(r_options_ptr, "maxit");
  if(r_opt != R_NilValue)
    opt.maxiter=(int)(REAL(r_opt)[0]);

  r_opt = getListElement(r_options_ptr, "size");
  if(r_opt != R_NilValue)
    opt.s=(int)(REAL(r_opt)[0]);

  r_opt = getListElement(r_options_ptr, "iprint");
  if(r_opt != R_NilValue)
    opt.IPrint=(int)(REAL(r_opt)[0]);

  r_opt = getListElement(r_options_ptr, "social");
  if(r_opt != R_NilValue)
    opt.nu=REAL(r_opt)[0];

  r_opt = getListElement(r_options_ptr, "tol");
  if(r_opt != R_NilValue)
    opt.tol=REAL(r_opt)[0];

  r_opt = getListElement(r_options_ptr, "delta");
  if(r_opt != R_NilValue)
    opt.delta=REAL(r_opt)[0];

  r_opt = getListElement(r_options_ptr, "ddelta");
  if(r_opt != R_NilValue)
    opt.ddelta=REAL(r_opt)[0];

  r_opt = getListElement(r_options_ptr, "idelta");
  if(r_opt != R_NilValue)
    opt.idelta=REAL(r_opt)[0];
  

  r_outf = getListElement(r_options_ptr, "outputfcn");
  opt.outfcn=&r_outfcn;

  r_opt = getListElement(r_options_ptr, "vectorized");
  if(r_opt != R_NilValue)
    opt.vectorized=(int)(REAL(r_opt)[0]);

  
  if (!setjmp(Jb)){
    signal(SIGFPE, catchfpe);
    
    exit_code=PSwarm(n, &r_objfun, lb, ub, lincons, A, b, &sol, &f, x0);
    

	/*if(opt.IPrint>=0){
	    if(exit_code){
			Rprintf("Abnormal exit\n");
		} else {
			Rprintf("Normal exit\n");
		}
	}*/	
	
  }

  /* Build a list for calling the objective function */
  PROTECT(s = allocVector(VECSXP, 3));
  PROTECT(t = allocVector(STRSXP, 3));
  
  r_fobj=createRRealScalar(f);
  r_x=createRRealVector(n,sol);
  r_ret=createRIntScalar(exit_code);
  
  SET_VECTOR_ELT(s, 0, r_x);
  SET_VECTOR_ELT(s, 1, r_fobj);
  SET_VECTOR_ELT(s, 2, r_ret);

  SET_VECTOR_ELT(t, 0, mkChar("x"));
  SET_VECTOR_ELT(t, 1, mkChar("f"));
  SET_VECTOR_ELT(t, 2, mkChar("ret"));

  setAttrib(s, R_NamesSymbol, t);
  
  UNPROTECT(2);

  if(lb)
    free(lb);
  if(ub)
    free(ub);
  if(A)
    free(A);
  if(b)
    free(b);
  if(x0)
    free(x0);
  if(sol)
    free(sol);
  
  return s;
  
}


