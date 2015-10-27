*deck @(#)mcvec.f	1.1  11/30/90
      subroutine mcvec(nsym,nbf,nob,nocc,cv,locsym,lok,len,mix,del,
     $     id2h,cr,nwwp, nfo)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcvec.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
      integer wpadti
      dimension nbf(2),nob(2),nocc(2),cv(2),nfo(2)
      dimension locsym(2),lok(2),len(2),mix(2),del(2)
      dimension cr(2)
c
c...
c...  real*8 cr(nwwp)
c...
      common /io/ inp,iout
c
c
      nob2 = 0
      do 100 l = 1, nsym
 100  nob2 = max(nob2,nob(l)**2)
c
      lhold = 1
      lu = lhold + nob2
      last = wpadti(lu + nob2)
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(last,cr,ncore,'mc vectors',0)
c
cpsdyn      if (last .lt. ncore) go to 190
cpsdyn      write (iout,9000)
cpsdyn 9000 format(' need more core for orbital transformation  bug!')
cpsdyn      call lnkerr(' ') 1234
c
 190  lla = 1
      llc = 1
c
      ll=0
c
      do 200 l = 1, nsym
         lla = lla + nfo(l)*nbf(l)
         if (nob(l) .eq. 0) go to 200
         llb = locsym(l) + 1
         call mcfrmu(cr(lhold),cr(lu),del(llc),lok(llb),len(llb),mix,
     $        nob(l),nocc(l),ll)
         call mcbcku(cr(lhold),cr(lu),cv(lla),nob(l),nbf(l),l,id2h)
         lla = lla + nob(l) * nbf(l)
         llc = llc + ll
 200  continue
c
      return
      end
