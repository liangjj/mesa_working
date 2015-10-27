#RUN = -execute -ce
RUN =

#######################################################################
#  This makefile runs the test programs for the linear equation routines
#  and the eigenvalue routines in LAPACK.  The test output files
#  are grouped as follows:
#
#       SLINTST,SEIGTST  -- Single precision real test routines
#       CLINTST,CEIGTST  -- Single precision complex test routines
#       DLINTST,DEIGTST  -- Double precision real test routines
#       ZLINTST,ZEIGTST  -- Double precision complex test routines
#
#  Test programs can be executed for all or some of the four different
#  precisions.  First, modify the RUN definition to refer to the
#  execute and desired execute options for your machine.
#  Then enter 'make' followed by one or more of the data types desired.
#  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments runs all eight test programs.
#  The executable files are called:
#       xlintsts, xlintstd, xlintstc, and xlintstz for LIN
#       xeigtsts, xeigtstd, xeigtstc, and xeigtstz for EIG
#  and exist in the current directory level.
#
#  To remove the output files after the tests have been run, enter
#       make clean
#
#  To re-run specific tests after a make, enter (for example):
#       'rm ssvd.out; make'  or:
#       'make ssvd.out' or:
#       'touch svd.in; make' (to re-run the single precision SVD tests.)
#
#       'rm *svd.out; make'  (to re-run all the SVD tests.)
#
#----------------------------------------------------------------------
#
#  Last updated by Susan Ostrouchov: March 7, 1993
#
#######################################################################

all:     single complex double complex16

SEIGTST= snep.out \
         ssep.out \
         ssvd.out \
         sec.out \
         sed.out \
         sgg.out \
         ssb.out \
	 ssg.out \
         sbal.out \
         sbak.out \
         sglm.out \
         sgqr.out \
         sgsv.out \
         slse.out

CEIGTST= cnep.out \
         csep.out \
         csvd.out \
         cec.out \
         ced.out \
         cgg.out \
         csb.out \
	 csg.out \
         cbal.out \
         cbak.out \
         cglm.out \
         cgqr.out \
         cgsv.out \
         clse.out

DEIGTST= dnep.out \
         dsep.out \
         dsvd.out \
         dec.out \
         ded.out \
         dgg.out \
         dsb.out \
	 dsg.out \
         dbal.out \
         dbak.out \
         dglm.out \
         dgqr.out \
         dgsv.out \
         dlse.out

ZEIGTST= znep.out \
         zsep.out \
         zsvd.out \
         zec.out \
         zed.out \
         zgg.out \
         zsb.out \
	 zsg.out \
         zbal.out \
         zbak.out \
         zglm.out \
         zgqr.out \
         zgsv.out \
         zlse.out


SLINTST= stest.out

CLINTST= ctest.out

DLINTST= dtest.out

ZLINTST= ztest.out

single:         $(SLINTST) $(SEIGTST)
complex:        $(CLINTST) $(CEIGTST)
double:         $(DLINTST) $(DEIGTST)
complex16:      $(ZLINTST) $(ZEIGTST)

#
# ======== SINGLE LIN TESTS ===========================

stest.out: stest.in xlintsts
	@echo Data file for testing REAL LAPACK linear equation routines
	$(RUN) xlintsts < stest.in > $@ 2>&1
#
# ======== COMPLEX LIN TESTS ==========================

ctest.out: ctest.in xlintstc
	@echo Data file for testing COMPLEX LAPACK linear equation routines
	$(RUN) xlintstc < ctest.in > $@ 2>&1
#
# ======== DOUBLE LIN TESTS ===========================

dtest.out: dtest.in xlintstd
	@echo Data file for testing DOUBLE PRECISION LAPACK linear equation routines
	$(RUN) xlintstd < dtest.in > $@ 2>&1
#
# ======== COMPLEX16 LIN TESTS ========================

ztest.out: ztest.in xlintstz
	@echo Data file for testing COMPLEX16 LAPACK linear equation routines
	$(RUN) xlintstz < ztest.in > $@ 2>&1
#
#
# ======== SINGLE EIG TESTS ===========================
#

snep.out: nep.in xeigtsts
	@echo NEP: Data file for testing Nonsymmetric Eigenvalue Problem routines
	$(RUN) xeigtsts < nep.in > $@ 2>&1

ssep.out: sep.in xeigtsts
	@echo SEP: Data file for testing Symmetric Eigenvalue Problem routines
	$(RUN) xeigtsts < sep.in > $@ 2>&1

ssvd.out: svd.in xeigtsts
	@echo SVD: Data file for testing Singular Value Decomposition routines
	$(RUN) xeigtsts < svd.in > $@ 2>&1

sec.out: sec.in xeigtsts
	@echo SEC: Data file for testing REAL Eigen Condition Routines
	$(RUN) xeigtsts < sec.in > $@ 2>&1

sed.out: sed.in xeigtsts
	@echo SEV: Data file for testing REAL Nonsymmetric Eigenvalue Driver
	$(RUN) xeigtsts < sed.in > $@ 2>&1

sgg.out: sgg.in xeigtsts
	@echo SGG: Data file for testing REAL Nonsymmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtsts < sgg.in > $@ 2>&1

ssb.out: ssb.in xeigtsts
	@echo SSB: Data file for testing REAL Symmetric Eigenvalue Problem routines
	$(RUN) xeigtsts < ssb.in > $@ 2>&1

ssg.out: ssg.in xeigtsts
	@echo SSG: Data file for testing REAL Symmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtsts < ssg.in > $@ 2>&1

sbal.out: sbal.in xeigtsts
	@echo SGEBAL: Data file for testing the balancing of a REAL general matrix
	$(RUN) xeigtsts < sbal.in > $@ 2>&1

sbak.out: sbak.in xeigtsts
	@echo SGEBAK: Data file for testing the back transformation of a REAL balanced matrix
	$(RUN) xeigtsts < sbak.in > $@ 2>&1

sglm.out: glm.in xeigtsts
	@echo GLM: Data file for testing Generalized Linear Regression Model routines
	$(RUN) xeigtsts < glm.in > $@ 2>&1

sgqr.out: gqr.in xeigtsts
	@echo GQR: Data file for testing Generalized QR and RQ factorization routines
	$(RUN) xeigtsts < gqr.in > $@ 2>&1

sgsv.out: gsv.in xeigtsts
	@echo GSV: Data file for testing Generalized Singular Value Decomposition routines
	$(RUN) xeigtsts < gsv.in > $@ 2>&1

slse.out: lse.in xeigtsts
	@echo LSE: Data file for testing Constrained Linear Least Squares routines
	$(RUN) xeigtsts < lse.in > $@ 2>&1
#
# ======== COMPLEX EIG TESTS ===========================

cnep.out: nep.in xeigtstc
	@echo NEP: Data file for testing Nonsymmetric Eigenvalue Problem routines
	$(RUN) xeigtstc < nep.in > $@ 2>&1

csep.out: sep.in xeigtstc
	@echo SEP: Data file for testing Symmetric Eigenvalue Problem routines
	$(RUN) xeigtstc < sep.in > $@ 2>&1

csvd.out: svd.in xeigtstc
	@echo SVD: Data file for testing Singular Value Decomposition routines
	$(RUN) xeigtstc < svd.in > $@ 2>&1

cec.out: cec.in xeigtstc
	@echo CEC: Data file for testing COMPLEX Eigen Condition Routines
	$(RUN) xeigtstc < cec.in > $@ 2>&1

ced.out: ced.in xeigtstc
	@echo CES: Data for the testing COMPLEX Nonsymmetric Schur Form Driver
	$(RUN) xeigtstc < ced.in > $@ 2>&1

cgg.out: cgg.in xeigtstc
	@echo CGG: Data file for testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtstc < cgg.in > $@ 2>&1

csb.out: csb.in xeigtstc
	@echo CHB: Data file for testing Hermitian Eigenvalue Problem routines
	$(RUN) xeigtstc < csb.in > $@ 2>&1

csg.out: csg.in xeigtstc
	@echo CSG: Data file for testing Symmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtstc < csg.in > $@ 2>&1

cbal.out: cbal.in xeigtstc
	@echo CGEBAL: Data file for testing the balancing of a COMPLEX general matrix
	$(RUN) xeigtstc < cbal.in > $@ 2>&1

cbak.out: cbak.in xeigtstc
	@echo CGEBAK: Data file for testing the back transformation of a COMPLEX balanced matrix
	$(RUN) xeigtstc < cbak.in > $@ 2>&1

cglm.out: glm.in xeigtstc
	@echo GLM: Data file for testing Generalized Linear Regression Model routines
	$(RUN) xeigtstc < glm.in > $@ 2>&1

cgqr.out: gqr.in xeigtstc
	@echo GQR: Data file for testing Generalized QR and RQ factorization routines
	$(RUN) xeigtstc < gqr.in > $@ 2>&1

cgsv.out: gsv.in xeigtstc
	@echo GSV: Data file for testing Generalized Singular Value Decomposition routines
	$(RUN) xeigtstc < gsv.in > $@ 2>&1

clse.out: lse.in xeigtstc
	@echo LSE: Data file for testing Constrained Linear Least Squares routines
	$(RUN) xeigtstc < lse.in > $@ 2>&1
#
# ======== DOUBLE EIG TESTS ===========================

dnep.out: nep.in xeigtstd
	@echo NEP: Data file for testing Nonsymmetric Eigenvalue Problem routines
	$(RUN) xeigtstd < nep.in > $@ 2>&1

dsep.out: sep.in xeigtstd
	@echo SEP: Data file for testing Symmetric Eigenvalue Problem routines
	$(RUN) xeigtstd < sep.in > $@ 2>&1

dsvd.out: svd.in xeigtstd
	@echo SVD: Data file for testing Singular Value Decomposition routines
	$(RUN) xeigtstd < svd.in > $@ 2>&1

dec.out: dec.in xeigtstd
	@echo DEC: Data file for testing DOUBLE PRECISION Eigen Condition Routines
	$(RUN) xeigtstd < dec.in > $@ 2>&1

ded.out: ded.in xeigtstd
	@echo DEV: Data file for DOUBLE PRECISION Nonsymmetric Eigenvalue Driver
	$(RUN) xeigtstd < ded.in > $@ 2>&1

dgg.out: dgg.in xeigtstd
	@echo DGG: Data file for testing DOUBLE PRECISION Nonsymmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtstd < dgg.in > $@ 2>&1

dsb.out: dsb.in xeigtstd
	@echo DSB: Data file for testing DOUBLE PRECISION Symmetric Eigenvalue Problem routines
	$(RUN) xeigtstd < dsb.in > $@ 2>&1

dsg.out: dsg.in xeigtstd
	@echo DSG: Data file for testing DOUBLE PRECISION Symmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtstd < dsg.in > $@ 2>&1

dbal.out: dbal.in xeigtstd
	@echo DGEBAL: Data file for testing the balancing of a DOUBLE PRECISION general matrix
	$(RUN) xeigtstd < dbal.in > $@ 2>&1

dbak.out: dbak.in xeigtstd
	@echo DGEBAK:  Data file for testing the back transformation of a DOUBLE PRECISION balanced matrix
	$(RUN) xeigtstd < dbak.in > $@ 2>&1

dglm.out: glm.in xeigtstd
	@echo GLM: Data file for testing Generalized Linear Regression Model routines
	$(RUN) xeigtstd < glm.in > $@ 2>&1

dgqr.out: gqr.in xeigtstd
	@echo GQR: Data file for testing Generalized QR and RQ factorization routines
	$(RUN) xeigtstd < gqr.in > $@ 2>&1

dgsv.out: gsv.in xeigtstd
	@echo GSV: Data file for testing Generalized Singular Value Decomposition routines
	$(RUN) xeigtstd < gsv.in > $@ 2>&1

dlse.out: lse.in xeigtstd
	@echo LSE: Data file for testing Constrained Linear Least Squares routines
	$(RUN) xeigtstd < lse.in > $@ 2>&1
#
# ======== COMPLEX16 EIG TESTS ===========================

znep.out: nep.in xeigtstz
	@echo NEP: Data file for testing Nonsymmetric Eigenvalue Problem routines
	$(RUN) xeigtstz < nep.in > $@ 2>&1

zsep.out: sep.in xeigtstz
	@echo SEP: Data file for testing Symmetric Eigenvalue Problem routines
	$(RUN) xeigtstz < sep.in > $@ 2>&1

zsvd.out: svd.in xeigtstz
	@echo SVD: Data file for testing Singular Value Decomposition routines
	$(RUN) xeigtstz < svd.in > $@ 2>&1

zec.out: zec.in xeigtstz
	@echo ZEC: Data file for testing COMPLEX16 Eigen Condition Routines
	$(RUN) xeigtstz < zec.in > $@ 2>&1

zed.out: zed.in xeigtstz
	@echo ZES: Data file for testing COMPLEX16 Nonsymmetric Schur Form Driver
	$(RUN) xeigtstz < zed.in > $@ 2>&1

zgg.out: zgg.in xeigtstz
	@echo ZGG: Data file for testing COMPLEX16 Nonsymmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtstz < zgg.in > $@ 2>&1

zsb.out: zsb.in xeigtstz
	@echo ZHB: Data file for testing Hermitian Eigenvalue Problem routines
	$(RUN) xeigtstz < zsb.in > $@ 2>&1

zsg.out: zsg.in xeigtstz
	@echo ZSG: Data file for testing Symmetric Generalized Eigenvalue Problem routines
	$(RUN) xeigtstz < zsg.in > $@ 2>&1

zbal.out: zbal.in xeigtstz
	@echo ZGEBAL: Data file for testing the balancing of a COMPLEX16 general matrix
	$(RUN) xeigtstz < zbal.in > $@ 2>&1

zbak.out: zbak.in xeigtstz
	@echo ZGEBAK: Data file for testing the back transformation of a COMPLEX16 balanced matrix
	$(RUN) xeigtstz < zbak.in > $@ 2>&1

zglm.out: glm.in xeigtstz
	@echo GLM: Data file for testing Generalized Linear Regression Model routines
	$(RUN) xeigtstz < glm.in > $@ 2>&1

zgqr.out: gqr.in xeigtstz
	@echo GQR: Data file for testing Generalized QR and RQ factorization routines
	$(RUN) xeigtstz < gqr.in > $@ 2>&1

zgsv.out: gsv.in xeigtstz
	@echo GSV: Data file for testing Generalized Singular Value Decomposition routines
	$(RUN) xeigtstz < gsv.in > $@ 2>&1

zlse.out: lse.in xeigtstz
	@echo LSE: Data file for testing Constrained Linear Least Squares routines
	$(RUN) xeigtstz < lse.in > $@ 2>&1
# ==============================================================================

xlintsts:
	cd LIN ; $(MAKE) single

xlintstc:
	cd LIN ; $(MAKE) complex

xlintstd:
	cd LIN ; $(MAKE) double

xlintstz:
	cd LIN ; $(MAKE) complex16

xeigtsts:
	cd EIG ; $(MAKE) single

xeigtstc:
	cd EIG ; $(MAKE) complex

xeigtstd:
	cd EIG ; $(MAKE) double

xeigtstz:
	cd EIG ; $(MAKE) complex16

clean:
	rm -f *.out core

cleanup:
	rm -f x* *.out core
