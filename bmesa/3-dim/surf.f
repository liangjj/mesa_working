*deck surf.f
c***begin prologue     surf
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           three-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            surface functions
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       surf
      subroutine surf(v3d,v2d,p,prj,srf,ind,nply,n3,nv,n2,nc,q,prnprj,
     1                prn2v,prn3v)
      implicit integer (a-z)
      real*8 v3d, v2d, prj, p, srf
      character*80 title
      character*2 q
      logical prnprj, prn2v, prn3v
      dimension v3d(n3,nv), v2d(n2,nc), p(nply)
      dimension prj(n2,nv), srf(nc,nv), ind(n3,4)
      common/io/inp, iout 
c     the character variable q can be q1, q2 or q3 depending on which
c     coordinate is normal to the surface.
c      if(prn2v) then
c         title='2d-eigenvectors'
c         call prntrm(title,v2d,n2,nc,n2,nc,iout)
c      endif
c      if(prn3v) then
c         title='3d-eigenvectors'
c         call prntrm(title,v3d,n3,nv,n3,nv,iout)
c      endif
      call rzero(prj,n2*nv)
      cnt3=0
      do 10 i=1,nply
         cnt2=0
         do 20 j=1,nply
            do 30 k=1,nply
               cnt2=cnt2+1
               cnt3=cnt3+1
               indx=ind(cnt3,4)
               do 40 l=1,nv  
                  prj(cnt2,l) = prj(cnt2,l) + v3d(indx,l)*p(i)
 40            continue
 30         continue
 20      continue
 10   continue
      if(prnprj) then
         title='primitive surface functions'   
         call prntrm(title,prj,n2,nv,n2,nv,iout)
      endif
      call ebtc(srf,v2d,prj,nc,n2,nv)
      if(prnprj) then
         title='surface functions'   
         call prntrm(title,srf,nc,nv,nc,nv,iout)
      endif                  
      return
      end       

