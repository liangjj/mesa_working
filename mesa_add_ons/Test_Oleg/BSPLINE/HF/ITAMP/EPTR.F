*
*     ------------------------------------------------------------------
*                       e p t r
*     ------------------------------------------------------------------
*
*       determines the position of the electron in the electron list
*
      SUBROUTINE eptr(el,elsymb, iel, *)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      CHARACTER el(*)*3, elsymb*3, bl*3
      data bl/'   '/
*
* ***** search electron list for lsymb
*
      if ( elsymb .eq. bl ) then
         iel = 0
         return
      endif
      do 10 i=1,nwf
         if (el(i) .eq. elsymb ) then
            iel = i
            return
         endif
10    continue
      iel = -1
      write(iscw,20) elsymb
20    format(/10x,a3,' not found in electron list')
      return 1
      end
