*deck @(#)fixbuf.f	5.1  11/6/94
      subroutine fixbuf(iwork,mrs,naoi,ncoi,icc,iac)
c
c***begin prologue     fixbuf
c***date written       871022   (yymmdd)
c***revision date      871116   (yymmdd)
c   16 november 1987   bhl      lanl
c   error checking and more sensible limits imposed on buffer space
c   correction are denoted c..bhl..limit
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)fixbuf.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       fixbuf
c
c
      implicit real*8 (a-h,o-z)
      if(ncoi.eq.0) then
c..bhl..limit         iac=iwork/(naoi*mrs)
         icc=0
         iac=0
c..bhl..limit         iac=iac*naoi
c..bhl..limit         iac=min(iac,naoi*ncoi)
      else
         iac=1
         ileft=iwork-naoi*mrs
         icc=ileft/mrs
         if(icc.eq.0) call lnkerr('  bug in mcscf fixbuf ')
         ic=iwork/(naoi*mrs)
         iac=ic/2
c..bhl..limit.start
         iac=min(iac,ncoi)
c..bhl..limit.end
         icc=(iwork-iac*naoi*mrs)/mrs
         iac=iac*naoi
c..bhl..limit         iac=min(iac,naoi*ncoi)
         icc=min(icc,ncoi*(ncoi+1)/2)
         if(icc.eq.0) then
            call lnkerr(' icc=0  bug in mcscf fixbuf..pass 2 ')
         end if
         if(iac.eq.0.and.naoi.ne.0) then
            call lnkerr(' iac=0  bug in mcscf fixbuf..pass 2 ')
         end if
      endif
c
c
      return
      end
