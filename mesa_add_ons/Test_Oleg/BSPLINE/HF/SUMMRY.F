*     ------------------------------------------------------------------
*               s u m m r y
*     ------------------------------------------------------------------
*
*       the results of a calculation are summarized.   these include
*   the following for each electron:
*
*          e(nl)   - diagonal energy parameter
*          i(nl)   - -(1/2)<nl|l|nl>
*          ke      - i(nl) + z <r>
*          rel     - relativistic shift (mass-velocity, darwin term,
*                    spin-spin contact term)
*          sigma   - screening parameter as defined by eq. (6-  ).
*          az(nl)  - starting parameter, p(r)/r**(l+1) as r -> 0.
*          1/r**3  - expected value of <1/r**3>
*          1/r     - expected value of <1/r>
*          r       - expected mean radius
*          r**2    - expected value of <r**2>
*
*   these results are followed by:
*
*          kinetic energy (ek)
*          potential energy (ep) = et - en
*          ratio                 = ep/en
*          non- relativistic energy (et - erel)
*          relativistic shift (erel) for the state
*          total energy (et)
*
      SUBROUTINE summry(et,erel)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /coeff/coef(100),ijptr(5,5)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
      DIMENSION r1(nwd),rm1(nwd),ek(nwd)
*
      pi = acos(-1.d0)
      write (3,9) atom,term
9     format(/// 24x,'atom ',a6,3x,'term ',a6//
     :   2x,'nl',7x,'e(nl)',
     :   8x,'i(nl)',7x,'ke(nl)',8x,'rel(nl)',3x,'s(nl)',7x,'az(nl)')
      en = 0.d0
      rel = .false.
*
*  *****  compute and print one-electron parameters
*
      do 10 i = 1,nwf
        r1(i) = quadr(i,i,1)
        ek(i) = -0.5d0*hl(el,i,i,rel)
        rm1(i) = quadr(i,i,-1)
        ekinp = ek(i) + z*rm1(i)
        en = en+ sum(i)*ekinp
	print *, i,sum(i),ek(i),rm1(i),ekinp
        rh = 3*n(i)*n(i) - l(i)*(l(i) + 1)
        sc = z - 0.5d0*rh/r1(i)
        s(i) = sc
        rels = rlshft(i,i)
        azz = az(z,h,ks,l(i)+1,p(l(i)+2,i))
        write (3,15)el(i),e(i,i),ek(i),ekinp,rels,s(i),azz
15      format(1x,a3,f14.7,3f13.6,f8.3,f14.6)
10    continue
*
*  *****  compute moments
*
      write(3,8) 'delta(r)'
 8    format(//2x,'nl',6x,a8,5x,'1/r**3',7x,'1/r',9x,'r',8x,'r**2')
      do 11 i = 1,nwf
        rm3 = 0
        if (l(i) .ne. 0) rm3 = quadr(i,i,-3)
        rp2 = quadr(i,i,2)
        rz = 0.
        if ( l(i) .eq. 0) rz = az(z,h,ks,l(i)+1,p(l(i)+2,i))**2/(4.*pi)
        write(3,16) el(i),rz,rm3,rm1(i),r1(i),rp2
16      format(1x,a3,f14.3,f13.4,f11.5,f10.5,f11.5)
11    continue
31    etn = et - erel
      epotl = etn - en
      ratio = epotl/en
      write(iscw,26) etn,en,erel,epotl,et,ratio
      write(3,26) etn,en,erel,epotl,et,ratio
26    format(//5x,'total energy (a.u.)'/5x,'----- ------'/
     : 10x,' non-relativistic   ',f15.8,t50,'kinetic   ',f15.8/
     : 10x,' relativistic shift ',f15.8,t50,'potential ',f15.8/
     : 10x,' relativistic       ',f15.8,t50,'ratio     ',f15.9)
13    return
      end
