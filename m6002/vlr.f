*deck @(#)vlr.f	1.1 9/7/91
c***begin prologue     vlr
c***date written       901008   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           vlr, link 6002
c***author             schneider, barry (lanl)
c***source             m6002
c***purpose            calculate long range potential
c***references         none
c
c***routines called
c***end prologue       vlr
      subroutine vlr (vints,ylm,mom,nomom,index,grid,npnts,nstri,
     1                nptmx,mxmom)
      implicit integer (a-z)
      real *8 vints, ylm, mom, grid, rval
      dimension vints(nptmx,nstri), mom(mxmom,nstri)
      dimension grid(4,nptmx), ylm(npnts), nomom(nstri)
      dimension index(3,mxmom,nstri)
      do 10 i=1,nstri
         do 20 mmom=1,nomom(i)
            l=index(1,mmom,i)
            m=index(2,mmom,i)
            power=index(3,mmom,i)
            if ( mom(mmom,i).ne.0.d+00) then
                 call lgndre(ylm,grid,npnts,l,m)
                 do 30 pt=1,npnts
                    rval=grid(1,pt)*grid(1,pt)
     1                   +grid(2,pt)*grid(2,pt)
     2                   +grid(3,pt)*grid(3,pt)
                    rval=sqrt(rval)
                    rval=1.d+00/(rval**power)
                    vints(pt,i)=vints(pt,i)+mom(mmom,i)*ylm(pt)
     1                                                 *rval
   30            continue 
            endif                                     
   20    continue
   10 continue
      return
      end
