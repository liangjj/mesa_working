*deck @(#)sinv.f	5.1  11/6/94
      subroutine sinv(s,smhalf,u,eigval,t1,t2,num,nnp,triang,
     $                iprint)
c***begin prologue     sinv
c***date written       850601  (yymmdd)
c***revision date      860813  (yymmdd)
c        13 august 1986  pws at lanl
c             changed call of rsp to degrsp, which rotates degenerate
c             eigenvectors to a hopefully sensible angle.
c
c***keywords           matrix, invert, square root
c***author             saxe, paul (lanl)
c***source
c***purpose                                              -1/2
c                      vectorized matrix square root  v=w     .
c***description
c                      call sinv(s,smhalf,u,eigval,t1,t2,num,nnp,
c                                triang,iprint)
c
c                        s       input matrix (nnp).
c                        smhalf  output matrix (nnp)
c                        u       scratch (num,num).
c                        eigval  eigenvalues of s (num).
c                        t1      scratch (num,num).
c                        t2      scratch (num,num).
c                        num     matrix dimension.
c                        nnp     num*(num+1)/2
c                        triang  scratch (nnp).
c                        iprint  print flag:  2**0    print overlap matrix.
c                                             2**1    print eigenvectors of s.
c                                             2**2    print s**-1/2.
c
c***references
c***routines called    print(math),vmove(math),rsp(clams),
c                      vecout(math),vsqrt(math),zero(math),ebct(math),
c                      ebc(math),sqtotr(math),vinv(math)
c***end prologue       sinv
      implicit integer (a-z)
c
      real*8 s(nnp),smhalf(nnp),u(num,num),t1(num,num)
      real*8 t2(num,num),eigval(num),triang(nnp)
c
      common /io/     inp,iout
c
c     ----- print the overlap matrix if requested -----
c
      if (iprint.eq.1) then
         write (iout,101)
 101     format (//,' the a.o. overlap matrix ',/)
         call print (s,nnp,num,iout)
      end if
c
c     ----- diagonalize s and get eigenvalues and vectors -----
c
      call vmove(triang,s,nnp)
      call degrsp(num,nnp,triang,eigval,1,u,t1,t2)
c
c     ----- print eigenvectors and eigenvalues of s if requested -----
c
      if (iprint.eq.2) then
         write (iout,102)
 102     format (//,' eigenvectors and values of the overlap matrix',/)
         call vecout(u,eigval,num)
       end if
c
c     ----- form s**(-1/2) -----
c
      call vinv(eigval,eigval,num)
      call rzero(t1,num**2)
      do 11 i=1,num
         t1(i,i)=sqrt(eigval(i))
   11 continue
      call ebct(t2,t1,u,num,num,num)
      call ebc(t1,u,t2,num,num,num)
      call sqtotr(smhalf,t1,num,nnp)
c
c     ----- print s**(-1/2) if requested -----
c
      if (iprint.eq.4) then
         write (iout,103)
 103     format (//,' s**(-1/2) matrix ',/)
         call print (smhalf,nnp,num,iout)
      end if
c
c
      return
      end
