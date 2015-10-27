*deck @(#)fillfc.f	5.1  11/6/94
      subroutine fillfc(fc,l,natoms3,natoms,bvec,ff,cref,nvar,
     $                  maxpt,temp1,temp2,order,x,a,grad,
     $                  bgvec,t,lp,lm)
c***begin prologue     fillfc.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)fillfc.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fillfc.f
      implicit none
c     --- input variables -----
      integer natoms3,natoms,nvar,maxpt
c     --- input arrays (unmodified) ---
      real*8 cref(3,natoms)
c     --- input arrays (scratch) ---
      integer order(natoms3)
      real*8 l(natoms3,6),bvec(6,natoms3)
      real*8 ff(nvar,maxpt),x(6),a(6,natoms3)
      real*8 temp1(natoms3,natoms3),temp2(natoms3,natoms3)
      real*8 grad(natoms3),bgvec(6),t(natoms3,6)
      real*8 lp(natoms3,6),lm(natoms3,6)
c     --- output arrays ---
      real*8 fc(natoms3,natoms3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer atom,ix,iy,iz,nrt,irt,jrt,jjrt
      integer i,j,jj,ixyz,iixyz
      real*8 small,zero,one
      real*8 detx
      logical trouble,debug
c
      parameter(debug=.false.,zero=0.0d+00,one=1.0d+00,small=1.0d-10)
c
      common /io/ inp,iout
c
 1000 format(//,' l vectors before schmidt',/)
 1010 format(//,' l vectors after schmidt',/)
 1020 format(/' a matrix'/)
 1030 format(/' cartesian gradient ',/,(3f20.6))
 1040 format(//,' t vectors ixyz,order(ixyz)',2i8/)
 1050 format(/'  b vectors before leq1f ')
 1060 format(/'  x vectors after leq1f ')
 1070 format(/'  fc first part done ')
 1080 format(//,' t vectors ixyz,order(ixyz)',2i8/)
 1090 format(/'  b vectors before leq1f ')
 1100 format(/'  x vectors after leq1f ')
 1110 format(/'  fc second part done ')
c
c     --- make the l vectors which comprise the matrix in the
c         linear system of equations
      call rzero(l,natoms3*6)
      do 4 atom=1,natoms
         ix=(atom-1)*3+1
         iy=(atom-1)*3+2
         iz=(atom-1)*3+3
         l(ix,1)=one
         l(iy,2)=one
         l(iz,3)=one
    4 continue
      do 5 atom=1,natoms
         l(3*atom-1,4)= cref(3,atom)
         l(3*atom  ,4)=-cref(2,atom)
    5 continue
      do 6 atom=1,natoms
         l(3*atom-2,5)= cref(3,atom)
         l(3*atom  ,5)=-cref(1,atom)
    6 continue
      do 7 atom=1,natoms
         l(3*atom-2,6) =  cref(2,atom)
         l(3*atom-1,6) = -cref(1,atom)
    7 continue
c
c     --- fill up the order array which allows us to refer to the fixed
c         coordinates (order(1)-order(6)) or to the varied coordinates
c         (order(7)-order(3n))
      do 20 i=1,natoms3
         order(i)=i
   20 continue
      order(6)=8
      order(7)=6
      order(8)=7
      if(debug) then
         write(iout,1000)
         call matout(l,natoms3,6,natoms3,6,iout)
      endif
c
c     --- orthogonalize the l vectors.
      nrt=6
      call schmidt(l,temp1,temp2,nrt,natoms3)
      if(debug) then
         write(iout,1010)
         call matout(l,natoms3,6,natoms3,6,iout)
      endif
c
c     --- contract the l vectors to form the matrix which
c         comprises the left hand side of the linear systems
      do 30 j=1,nrt
         do 25 i=1,nrt
            a(j,i)=l(order(i),j)
   25    continue
   30 continue
      if(debug) then
         write(iout,1020)
         call matout(a,nrt,nrt,nrt,nrt,iout)
      endif
c
c     --- check for trouble.
      trouble=.false.
      do 40 j=1,nrt
         do 35 i=1,nrt
         if(fc(order(i),order(j)).gt.small)
     $                        trouble=.true.
   35    continue
   40 continue
c
c     --- solve for the grad vector
      do 50 jrt=1,nrt
         bgvec(jrt)=zero
         do 45 jj=nrt+1,natoms3
            bgvec(jrt)=bgvec(jrt) -
     $                 l(order(jj),jrt)*grad(order(jj))
   45    continue
   50 continue
      call scopy(nrt,bgvec,1,a(1,nrt+1),1)
      call flin(a,nrt,nrt,1,detx)
      call scopy(nrt,a(1,nrt+1),1,bgvec,1)
      do 60 i=1,nrt
         grad(order(i))=bgvec(i)
   60 continue
      if(debug) then
         write(iout,1030) (grad(i),i=1,natoms3)
      endif
c
c     --- loop over columns of the force constant matrix
c         form the b-vector and solve the linear system of
c         equations. it is necessary to loop first over the
c         non fixed coordinates since these force constant
c         matrix elements are needed to finish
c
c         ixyz runs over columns of the force constant matrix 
c         non-fixed coordinates                              
      do 100 ixyz=nrt+1,natoms3
c        --- irt runs over translations and rotations ---
         do 70 irt=1,nrt
            bvec(irt,ixyz)=zero
            do 65 j=nrt+1,natoms3
               bvec(irt,ixyz)=bvec(irt,ixyz) - l(order(j),irt) *
     $                                   fc(order(j),order(ixyz))
   65       continue
   70    continue
c
c     --- this is where the gradient dependent term is added
c         to the b vector
         call tvec(t,cref,order(ixyz),lp,lm,natoms3,natoms,
     $             temp1,temp2,nrt)
         if(debug) then
            write(iout,1040) ixyz,order(ixyz)
            call matout(t,natoms3,6,natoms3,6,iout)
         endif
         do 80 jjrt=1,nrt
            do 75 iixyz=1,natoms3
               bvec(jjrt,ixyz)=bvec(jjrt,ixyz) -
     $                        t(iixyz,jjrt)*grad(iixyz)
   75       continue
   80    continue
  100 continue
c
c     --- solve the linear equations for the missing elements of the 
c         force constant matrix                                      
      if(debug) then
         write(iout,1050)
         call matout(bvec,6,natoms3,6,natoms3,iout)
      endif
      do 150 j=1,nrt
         do 145 i=1,nrt
            a(j,i)=l(order(i),j)
  145    continue
  150 continue
      call scopy(nrt*natoms3,bvec,1,a(1,nrt+1),1)
      call flin(a,nrt,nrt,natoms3,detx)
      call scopy(nrt*natoms3,a(1,nrt+1),1,bvec,1)
      if(debug) then
         write(iout,1060)
         call matout(bvec,6,natoms3,6,natoms3,iout)
      endif
      do 160 ixyz=nrt+1,natoms3
         do 155 i=1,nrt
            fc(order(i),order(ixyz))=bvec(i,ixyz)
            fc(order(ixyz),order(i))=bvec(i,ixyz)
  155    continue
  160 continue
      if(debug) then
         write(iout,1070)
         call matout(fc,natoms3,natoms3,natoms3,natoms3,iout)
      endif
c
c     --- ixyz runs over columns of the force constant matrix 
c         fixed coordinates                                   
      do 200 ixyz=1,nrt
c        --- irt runs over translations and rotations 
         do 170 irt=1,nrt
            bvec(irt,ixyz)=zero
            do 165 j=nrt+1,natoms3
               bvec(irt,ixyz)=bvec(irt,ixyz) - l(order(j),irt) *
     $                                    fc(order(j),order(ixyz))
  165       continue
  170    continue
c
c        --- this is where the gradient dependent term is added 
c            to the b vector
         call tvec(t,cref,order(ixyz),lp,lm,natoms3,natoms,
     $             temp1,temp2,nrt)
         if(debug) then
            write(iout,1080) ixyz,order(ixyz)
            call matout(t,natoms3,6,natoms3,6,iout)
         endif
         do 180 jjrt=1,nrt
            do 175 iixyz=1,natoms3
               bvec(jjrt,ixyz)=bvec(jjrt,ixyz) -
     $                        t(iixyz,jjrt)*grad(iixyz)
  175       continue
  180    continue
  200 continue
c
c     --- solve the linear equations for the missing elements of the 
c     force constant matrix                                      
      do 250 j=1,nrt
         do 245 i=1,nrt
            a(j,i)=l(order(i),j)
  245    continue
  250 continue
      if(debug) then
         write(iout,1090)
         call matout(bvec,6,natoms3,6,natoms3,iout)
      endif
      call scopy(nrt*nrt,bvec,1,a(1,nrt+1),1)
      call flin(a,nrt,nrt,nrt,detx)
      call scopy(nrt*nrt,a(1,nrt+1),1,bvec,1)
      if(debug) then
         write(iout,1100)
         call matout(bvec,6,natoms3,6,natoms3,iout)
      endif
      do 260 ixyz=1,nrt
         do 255 i=1,nrt
            fc(order(i),order(ixyz))=bvec(i,ixyz)
            fc(order(ixyz),order(i))=bvec(i,ixyz)
  255    continue
  260 continue
      if(debug) then
         write(iout,1110)
         call matout(fc,natoms3,natoms3,natoms3,natoms3,iout)
      endif
c
c
      return
      end
