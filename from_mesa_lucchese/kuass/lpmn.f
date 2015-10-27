        SUBROUTINE LPMN(MM,M,N,PM,narg,xarg)
C        SUBROUTINE LPMN(MM,M,N,X,PM,PD,narg,xarg)
C
C       =====================================================
C       Purpose: Compute the associated Legendre functions 
C                Pmn(x) and their derivatives Pmn'(x)
C       Input :  x  --- Argument of Pmn(x)
C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
C                mm --- Physical dimension of PM and PD
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:MM,0:MM,*),xarg(*)
c        DIMENSION PM(0:MM,0:N),PD(0:MM,0:N)
        DO 100 IARG=1,NARG
           X=XARG(IARG)
        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I,IARG)=0.0D0
c10         PD(J,I,IARG)=0.0D0
 10        continue
        PM(0,0,IARG)=1.D0   
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I,IARG)=X**I
 15           CONTINUE
        else
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I,IARG)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1,IARG)
        DO 35 I=0,M
35         PM(I,I+1,IARG)=(2.0D0*I+1.0D0)*X*PM(I,I,IARG)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J,iarg)=((2.0D0*J-1.0D0)*X*PM(I,J-1,IARG)-
     &             (I+J-1.0D0)*PM(I,J-2,IARG))/(J-I)
40      CONTINUE
        endif
 100    CONTINUE
        RETURN
        END
