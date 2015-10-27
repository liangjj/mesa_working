*deck @(#)schmdt.f	1.1 9/8/91
c***begin prologue     schmdt
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gram schmidt
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            schmidt orthogonalize a set of orbitals. 
c***
c***                 
c***                 
c***references       
c
c***routines called    cebtc(mylib),cvmul(math),iosys(io),sscal(clams)
c***end prologue       m7001
      subroutine schmdt(fns,ddfns,fnsc,ddfnsc,mask,ovlpr,ovlpc,work,n,
     1                  nbf,type)
      implicit integer (a-z)
      real *8 fns, ddfns, mask, sdot, ovlpr, tol, work
      complex *16 fnsc, ddfnsc, cdotu, ovlpc
      character *(*) type
      dimension fns(n,nbf), ddfns(n,nbf), ovlpr(nbf)
      dimension  fnsc(n,nbf), ovlpc(nbf), mask(nbf,nbf)
      dimension ddfnsc(n,nbf), work(nbf)
      parameter (tol=1.d-08)
      common /io/ inp, iout
c**********************************************************************c
c                   normalize the first function                       c
c**********************************************************************c
      if (type.eq.'real') then
          do 10 i=1,nbf
             work(i)=0.d0
             ovlpr(1)=sdot(n,fns(1,i),1,fns(1,i),1)
             if (ovlpr(1).gt.tol) then
                 work(i)=1.d0
                 ovlpr(1)=1.d0/sqrt(ovlpr(1))
                 call sscal(n,ovlpr(1),fns(1,i),1)
                 call sscal(n,ovlpr(1),ddfns(1,i),1)
                 if (i.gt.1) then
c**********************************************************************c
c              overlaps with previous functions                        c
c**********************************************************************c
                     call ebtc(ovlpr,fns,fns(1,i),i-1,n,1)
                     do 20 j=1,i-1
                        ovlpr(j)=ovlpr(j)*mask(i,j)
   20                continue         
c**********************************************************************c
c               subtract off overlaps                                  c
c**********************************************************************c
                     call ambc(fns(1,i),fns,ovlpr,n,i-1,1)
                     call ambc(ddfns(1,i),ddfns,ovlpr,n,i-1,1)
c**********************************************************************c
c                normalize                                             c
c**********************************************************************c
                     ovlpr(1)=sdot(n,fns(1,i),1,fns(1,i),1)
                     work(i)=0.d0
                     if (ovlpr(1).gt.tol) then
                         work(i)=1.d0 
                         ovlpr(i)=1.d0/sqrt(ovlpr(1))
                         call sscal(n,ovlpr(i),fns(1,i),1)
                         call sscal(n,ovlpr(i),ddfns(1,i),1)
                     endif
                 endif
             endif
   10     continue     
c**********************************************************************c
c                   normalize the first function                       c
c**********************************************************************c
      else
          do 30 i=1,nbf
             work(i)=0.d0
             ovlpc(1)=cdotu(n,fnsc(1,i),1,fnsc(1,i),1)
             if (abs(ovlpc(1)).gt.tol) then
                 work(i)=1.d0
                 ovlpc(1)=1.d0/sqrt(ovlpc(1))
                 call cscal(n,ovlpc(1),fnsc(1,i),1)
                 call cscal(n,ovlpc(1),ddfnsc(1,i),1)
                 if (i.gt.1) then
c**********************************************************************c
c              overlaps with previous functions                        c
c**********************************************************************c
                     call cebtc(ovlpc,fnsc,fnsc(1,i),i-1,n,1)
                     do 40 j=1,i-1
                        ovlpc(j)=ovlpc(j)*mask(i,j)
   40                continue         
c**********************************************************************c
c               subtract off overlaps                                  c
c**********************************************************************c
                     call cambc(fnsc(1,i),fnsc,ovlpc,n,i-1,1)
                     call cambc(ddfnsc(1,i),ddfnsc,ovlpc,n,i-1,1)
c**********************************************************************c
c                normalize                                             c
c**********************************************************************c
                     ovlpc(1)=cdotu(n,fnsc(1,i),1,fnsc(1,i),1)
                     work(i)=0.d0
                     if (abs(ovlpc(1)).gt.tol) then
                         work(i)=1.d0 
                         ovlpc(1)=1.d0/sqrt(ovlpc(1))
                         call cscal(n,ovlpc(1),fnsc(1,i),1)
                         call cscal(n,ovlpc(1),ddfnsc(1,i),1)
                     endif
                 endif
             endif
   30     continue     
      endif
      return
      end
