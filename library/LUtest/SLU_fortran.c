/*    
      Bridge from fortran to c for calling SuperLU routines
      Written by Mark Baertschy starting from the similar 
      routine provided with SuperLU.  Allows ordinary fortran
      indexing with first element array(1).

      NOTE: The principal feature of the bridges 
      is a choice of operations, in particular the ability 
      to save the factorization so that subsequent calls to
      solve the equations with different right hand sides
      run much faster.  For that reason the routine must be 
      called to release memory after its use.  Subsequent calls
      without doing this simply allocate more memory.

*/

#include <stdlib.h>
#include <stdio.h>

#include "zsp_defs.h"
#include "util.h"
#include "Cnames.h"


int
slu_fortran__(int *n, int *nnz, int *nrhs, doublecomplex *values,
                int *rowind, int *colptr, int *perm_c,
                doublecomplex *b, int *ldb, int *info, char *firsttime_in)

/* firsttime:  'Y'-->reorder, factorize, solve
            'N'-->just solve (reuse old factorization)
            'P'-->just factor (no solve)
            'F'-->free up memory (discard old factorization) */
{
    static SuperMatrix L, U;
    SuperMatrix A, B;
    static SCformat *Lstore;
    static NCformat *Ustore;
    static int      *perm_r; /* row permutations from partial pivoting */
    static int      panel_size, permc_spec, i;
/* variables for test write statements below
    static int  imccurdy, imcc;
*/
    static mem_usage_t   mem_usage;
    char trans[1], firsttime;
    double *utime, t1;
    extern SuperLUStat_t SuperLUStat;
    int Fortran,Factor,Solve;

    firsttime=*firsttime_in;
    *trans = 'N';
    utime = SuperLUStat.utime;


    switch (firsttime) {
    case 'Y': case 'y':
      firsttime = 'Y';
      printf("SuperLU: factor and solve\n");
      break;
    case 'N': case 'n':
      firsttime = 'N';
      printf("SuperLU: solve, use old factorization\n");
      break;
    case 'P': case 'p':
      firsttime = 'P';
      printf("SuperLU: factor but no solve\n");
      break;
    case 'F': case 'f':
      firsttime = 'F';
      printf("SuperLU: delete old factorization\n"); 
      break;
    default:
      printf("SuperLU: unknown option for refact\n");
      return (-1);
    }


    Factor = ( ( (firsttime=='Y') || (firsttime=='P') ) && (firsttime!='F') );
    Solve = ( ( (firsttime=='Y') || (firsttime=='N') ) && (firsttime!='F') );
    Fortran = (colptr[0]==1);
    if (Factor&&Fortran) {
      /* Adjust to 0-based indexing */
      printf("SuperLU: adjusting indices\n");
      for (i = 0; i < *nnz; ++i) --rowind[i];
      for (i = 0; i <= *n; ++i) --colptr[i];
    }
    
    if (Factor) {
      zCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr, NC, _Z, GE);
    }
    if (Solve) {
      zCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, DN, _Z, GE);
    } else {
      zCreate_Dense_Matrix(&B, 1, 1, b, 1, DN, _Z, GE);
    }
    printf("Created B\n");

    if (Factor) {
      if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
    }
    if (Factor) {
      /*
       * Get column permutation vector perm_c[], according to permc_spec:
       *   permc_spec = 0: use the natural ordering 
       *   permc_spec = 1: use minimum degree ordering on structure of A'*A
       *   permc_spec = 2: use minimum degree ordering on structure of A'+A
       */
      permc_spec = 2;
      printf("SuperLU: getting permutation.\n");
      get_perm_c(permc_spec, &A, perm_c);
    }
    *info = 0;
    
    if (Factor) {
      /*zgssv(&A, perm_c, perm_r, &L, &U, &B, info, Solve);*/
      int relax;
      char refact[1];
      int lwork = 0, *etree;
      SuperMatrix AC;
      double diag_pivot_thresh = 0.0;
      double drop_tol = 0;
      DNformat *Bstore;

      *refact = 'N';
      panel_size = sp_ienv(1);
      relax = sp_ienv(2);
      Bstore = B.Store;

      if ( A.nrow != A.ncol || A.nrow < 0 ||
           A.Stype != NC || A.Dtype != _Z || A.Mtype != GE )
          *info = -1;

      if ( *info != 0 ) {
          i = -(*info);
          xerbla_("zgssv", &i);
          return;
      }

      StatInit(panel_size, relax);
      utime = SuperLUStat.utime;

      if ( !(etree = intMalloc(A.ncol)) ) ABORT("Malloc fails for etree[].");

      t1 = SuperLU_timer_();
      sp_preorder(refact, &A, perm_c, etree, &AC);
      utime[ETREE] = SuperLU_timer_() - t1;

      t1 = SuperLU_timer_();
      zgstrf(refact, &AC, diag_pivot_thresh, drop_tol, relax, panel_size,
             etree, NULL, lwork, perm_r, perm_c, &L, &U, info);
      utime[FACT] = SuperLU_timer_() - t1;

      SUPERLU_FREE (etree);
      Destroy_CompCol_Permuted(&AC);

      PrintStat( &SuperLUStat );
      StatFree();

      if ( *info == 0 ) {
          panel_size = sp_ienv(1);

	  Lstore = (SCformat *) L.Store;
	  Ustore = (NCformat *) U.Store;
	  printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
	  printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
	  printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
	
	  zQuerySpace(&L, &U, panel_size, &mem_usage);
	  printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		 mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
		 mem_usage.expansions);
	
      } else {
          printf("zgssv() error returns INFO= %d\n", *info);
	  if ( info <= n ) { /* factorization completes */
	    zQuerySpace(&L, &U, panel_size, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
		   mem_usage.expansions);
	  }
      }
      Destroy_SuperMatrix_Store(&A);
      printf("Done with factorization\n");
    } 
    if ((*info == 0) && Solve) {
        if ( B.ncol < 0 ) {
            *info = -6;
            printf("problem with B.ncol\n");
        }
        if ( B.Stype != DN ) {
            *info = -6;
            printf("problem with B.Stype\n");
        }
        if ( B.Dtype != _Z ) {
            *info = -6;
            printf("problem with B.Dtype\n");
        }
        if ( B.Mtype != GE ) {
            *info = -6;
            printf("problem with B.Mtype\n");
        }
        t1 = SuperLU_timer_();
	/* Solve the system A*X=B, overwriting B with X. */
    printf("Proceeding to solve\n");
	zgstrs (trans, &L, &U, perm_r, perm_c, &B, info);
    printf("done with solve, info = %i\n",*info);
        utime[SOLVE] = SuperLU_timer_() - t1;
    }

    if ( firsttime == 'F' ) {
      SUPERLU_FREE (perm_r);
      Destroy_SuperNode_Matrix(&L);
      Destroy_CompCol_Matrix(&U);
    }

    Destroy_SuperMatrix_Store(&B);
    if (Factor&&Fortran) {
      /* Restore to 1-based indexing */
      for (i = 0; i < *nnz; ++i) ++rowind[i];
      for (i = 0; i <= *n; ++i) ++colptr[i];
    }
    return 0;
}
