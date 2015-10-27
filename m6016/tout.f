      subroutine tout(t,tij,ki,i,j,nchan,ntchn,ni,nj,nlm)
      implicit integer (a-z)
      complex*16  t, tij
      character *3 itoc
      character *80 title
      real *8 cross, sigtot, ki, pi
      dimension t(ntchn,ntchn), tij(ni,nj), nlm(nchan)
      common /io/ inp,iout
      data pi / 3.141592653589793  /
      title='tij-matrix for i= '//itoc(i)//' j= '//itoc(j)
      iloc=0
      do 10 ii=1,i-1
         iloc=iloc+nlm(ii)
   10 continue
      jloc=0
      do 20 jj=1,j-1
         jloc=jloc+nlm(jj)
   20 continue
      icnt=iloc
      do 30 ii=1,ni
         icnt=icnt+1
         jcnt=jloc
         do 40 jj=1,nj
            jcnt=jcnt+1
            tij(ii,jj)=t(icnt,jcnt)
   40    continue
   30 continue
      cross=0.d+00
      sigtot=0.d+00
      do 50 ii=1,ni
         do 60 jj=1,nj
            cross=cross+abs(tij(ii,jj))**2
            sigtot=sigtot+imag(tij(ii,jj))
   60    continue
   50 continue
      cross=4.d+00*pi*cross/(ki*ki)  
      sigtot=4.d+00*pi*sigtot/(ki*ki)
      call prntcm(title,tij,ni,nj,ni,nj,iout) 
      write (iout,70) cross
      write (iout,80) sigtot
   70 format (//,5x,'total cross section from standard formula',
     1               1x,f20.10)
   80 format (//,5x,'total cross section from optical theorem ',
     1               1x,f20.10)
    1 format (/,5x,a80)
      return
      end
