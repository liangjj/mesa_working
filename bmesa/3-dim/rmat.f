*deck rmat.f
c***begin prologue     rmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           r-matrix
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form partial three-dimensional r-matrix from 
c***                   surface projections on two planes.  
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       rmat
      subroutine rmat(srfi,srfj,eig,energy,rmt,n2,n3,nen,ic,jc,prnt)
      implicit integer (a-z)
      real*8 srfi, srfj, eig, energy, rmt, half, ehalf
      character*80 title
      character*16 fptoc
      character*1 itoc
      logical prnt
      dimension srfi(n2,n3), srfj(n2,n3), eig(n3), energy(nen)
      dimension rmt(n2,n2)
      common/io/inp, iout 
      data half/.5d0/
c      title='i surface'
c      call prntrm(title,srfi,n2,n3,n2,n3,iout)
c      title='j surface'
c      call prntrm(title,srfj,n2,n3,n2,n3,iout)
      if(ic.eq.jc) then
         do 10 ene=1,nen
            call rzero(rmt,n2*n2)
            ehalf = half*energy(ene)
            do 20 i=1,n2
               do 30 j=1,i
                  do 40 k=1,n3
                     rmt(i,j) = rmt(i,j) + 
     1                             srfi(i,k)*srfj(j,k)/(eig(k)-ehalf)
 40               continue
                  rmt(i,j)=half*rmt(i,j)
                  rmt(j,i) = rmt(i,j)  
 30            continue
 20         continue
            if(prnt) then
               title='r-matrix for surface i = '//itoc(ic)//
     1               ' surface j = '// itoc(jc)//
     2               ' at energy = '//fptoc(ehalf)
               call prntrm(title,rmt,n2,n2,n2,n2,iout)
            endif
 10      continue   
      else
         do 50 ene=1,nen
            call rzero(rmt,n2*n2)
            ehalf = half*energy(ene)
            do 60 i=1,n2
               do 70 j=1,n2
                  do 80 k=1,n3
                     rmt(i,j) = rmt(i,j) + 
     1                             srfi(i,k)*srfj(j,k)/(eig(k)-ehalf)
 80               continue
                  rmt(i,j)=half*rmt(i,j)
 70            continue
 60         continue
            if(prnt) then
               title='r-matrix for surface i = '//itoc(ic)//
     1               ' surface j = '// itoc(jc)//
     2               ' at energy = '//fptoc(ehalf)
               call prntrm(title,rmt,n2,n2,n2,n2,iout)
            endif
 50      continue   
      endif
      return
      end       


