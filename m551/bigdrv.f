*deck @(#)bigdrv.f	1.1  11/30/90
      subroutine bigdrv(cr,icr,nwint,buf,
     $     thrv,thrn,thrs,
     $     n,nmax,nroot,nsml,nlook,msml,diag,lenb)
c
c***begin prologue     bigdrv
c***date written       871022   (yymmdd)
c***revision date      880130   (yymmdd)
c   30 january 1988    bhl at brl
c      definition of iblk changed to reflect the fact that nroot
c      new vectors are created each iteration
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)bigdrv.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       bigdrv
c
c
      implicit real*8(a-h,o-z)
      integer wpadti
c
c
      dimension cr(2),icr(2),buf(2),thrv(20)
c
      common /io/ inp,iout
c...
c...  equivalence (cr,icr)
c...  dimension icr(nwint)
c...
c
c
c     call printm(buf,5,1)
      mdim=max(nsml,nmax*msml+msml)
      itr=(mdim*(mdim+1))/2
      isqr=mdim*mdim
      iblk=nmax*n*nroot+msml*n
c
c     allocate core
      ib=1
      ic=ib+iblk
      ish=ic+iblk+msml*n
      idiag=ish+itr
      itemp=idiag+n
      ieval=itemp+itr
      ievec=ieval+mdim
      ipt=wpadti(ievec+isqr)
      iscr=iadtwp(ipt+nsml)
      need=wpadti(iscr+5*n)
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(need,icr,ngot,'davidson driver',0)
c
c
      call bigdav(cr(ic),cr(ib),cr(ish),
     $     diag,cr(itemp),cr(ieval),
     $     cr(ievec),buf,icr(ipt),
     $     cr(iscr),thrv,thrn,thrs,
     $     n,nmax,nroot,nsml,nlook,msml,lenb)
c
      return
      end
