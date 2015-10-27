*deck @(#)linvec.f	5.1  11/6/94
      subroutine linvec(s,smhalf,eigvec,eigval,t1,t2,numin,numout,
     1                  nnp,tol,dir,prnt,noprnt)
c***begin prologue     linvec
c***date written       910528  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           vectors, matrix, invert, square root
c***author             schneider, barry (lanl)
c***source
c***purpose                                              
c 
c***description
c                      call linvec(s,smhalf,u,eigval,t1,t2,numin,
c                                  numout,nnp,triang,tol,prnt,noprnt)
c
c                        s       input matrix (numin,numin)
c                        smhalf  output matrix (nnp) if requested
c                        eigvec  eigenvectors of s (numin,numin).
c                        eigval  eigenvalues of s (numin).
c                        t1      scratch (numin,numin).
c                        t2      scratch (numin,numin).
c                        numin   matrix dimension.
c                        numout number vectors above tol
c                        nnp     numin*(numin+1)/2
c                        tol     acceptance for eigenvalue
c                                ( 10**-6 reasonable )
c                        dir     if = s**-1/2 will compute same
c                        prnt    print flag:   overlap
c                                              eigenvectors
c                                              s**-1/2
c                        noprnt  number to print
c***references
c***routines called    print(math),rsp(clams),aeqbc(mylib)
c                      sqtotr(util)
c                      
c***end prologue       linvec
      implicit integer (a-z)
c
      real*8 s(numin,numin), smhalf(nnp), eigvec(numin,*)
      real*8 tmp
      real*8 t2(numin,numin), eigval(numin), t1(numin,numin), tol
      character*16 bflabl(200)
      character *(*) prnt, dir
c
      common /io/ inp,iout
c
c
      call locase(dir,dir)
c     ----- print the overlap matrix if requested -----
      call locase(prnt,prnt)      
c
      call sqtotr(smhalf,s,numin,nnp)
      if (prnt.eq.'overlap') then
          write (iout,101)
          call print (smhalf,nnp,numin,iout)
       end if
c
c     ----- diagonalize s(smhalf) and get eigenvalues and vectors -----
c
      call degrsp(numin,nnp,smhalf,eigval,1,eigvec,t1,t2)
      write(iout,102) (eigval(i),i=1,numin)
      if (prnt.eq.'eigenvectors') then
         call iosys('read character "basis function labels" from rwf',
     $              -1,0,0,bflabl)
         call wvec(eigvec,eigval,numin,noprnt,bflabl,' ')
      end if
c
c     ----- form s**(-1/2) -----
c     ----- using only eigenvectors above tol -----
      numout=0
      do 10 i=1,numin
         if ( eigval(i).lt.0.d+00 ) then
              write (iout,103) eigval(i)
         else if ( abs(eigval(i)).gt.tol ) then
              numout=numout+1
              eigval(numout)=eigval(i)
              call scopy(numin,eigvec(1,i),1,eigvec(1,numout),1)
         end if
   10 continue  
      write (iout,104) numout, numin
      call scopy(numout*numin,eigvec,1,t1,1)
c     ----- scale vectors by eigval to get orthonormal set -----
      do 20 i=1,numout
         tmp=1.e+00/sqrt(eigval(i)) 
         do 30 j=1,numin
            eigvec(j,i)=eigvec(j,i)*tmp
   30    continue
   20 continue
c     ----- eigvec now contains an orthornormal set of linearly -----
c     -----          independent vectors -----
c     form s**-1/2
      if (dir.eq.'s**-1/2') then
          call rzero(t2,numin*numin)
          call aeqbc(t2,1,numin,t1,1,numin,eigvec,numin,1,
     1               numin,numout,numin)
          call sqtotr(smhalf,t2,numin,nnp)
c
c     ----- print s**(-1/2) if requested -----
c
          if ( prnt.eq.'s**-1/2') then
               write (iout,105)
               call print (smhalf,nnp,numin,iout)
          endif
      endif
c
c
      return
  101 format (//,' the a.o. overlap matrix ',/)
  102 format(/,' eigenvalues of the overlap matrix ',
     $           (/,5x,4(d15.8,1x)))
  103 format (/,5x,'negative eigenvalue of overlap',1x,d15.8)
  104 format (/,5x,i5,1x,'linearly independent vectors out of',1x,i5)
  105 format (//,' s**(-1/2) matrix ',/)
      end
