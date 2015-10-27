      subroutine grnset (gr1,gr2,en,rk,charge,pt,rmtrad,ncst,jind,chnloc
     1 ,nsts,ntchn,nptmx,lplsmx,ncmax,ops)
*
      implicit integer(a-z)
      character *(*) ops
      logical logkey
      real *8gr1, gr2, en, rk, charge, pt, rmtrad, cjn, enc, ark
      real *8f, fp, g, gp
      common /io/ inp, iout
      dimension gr1(nptmx,ntchn), gr2(nptmx,ntchn)
      dimension chnloc(lplsmx,nsts), en(ncmax,nsts), rk(ncmax,nsts)
      dimension ncst(nsts), pt(nptmx)
      dimension jind(ncmax,nsts)
*
*
*          the greens functions are defined as solutions
*          of the equation;
*            g'' + 2*abs(e) + l*(l+1)/(r+r) -2*v = delta(r-r')
*                -
*          the potential has already been mutiplied by
*          two in the readv routines. for a free particle
*          this leads to the greens function being the negative
*          of the product of gr1 and gr2 as defined in the code.
      do 30 is=1,nsts
      nsch=ncst(is)
      do 20 i=1,nsch
      lx=jind(i,is)-1
      cjn=0.d+00
      ncl=chnloc(jind(i,is),is)
*
      enc=en(i,is)
      ark=rk(i,is)
      call grncal (rmtrad,lx,enc,charge,f,fp,g,gp)
      cjn=-gp/fp
*
      do 10 ir=1,nptmx
      call grncal (pt(ir),lx,enc,charge,f,fp,g,gp)
      gr1(ir,ncl)=f
      gr2(ir,ncl)=(g+cjn*f)/rk(i,is)
   10 continue
   20 continue
   30 continue
*
*
      if (logkey(ops,'print=lam=green-fn',.false.,' ')) then
      write (iout,40)
      call matprt (gr1,nptmx,ntchn,nptmx,ntchn,0,0,0,0,0,0,0)
      write (iout,50)
      call matprt (gr2,nptmx,ntchn,nptmx,ntchn,0,0,0,0,0,0,0)
      endif
      return
*
c
   40 format (/,5x,' regular solution')
   50 format (/,5x,' irregular solution')
      end
