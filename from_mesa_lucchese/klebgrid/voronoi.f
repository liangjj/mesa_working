      subroutine voronoi(miter,ncent,a,ngrid,grid,npmx,d,r,p,wtt
     $,irad,rad,arad)
      implicit real*8 (a-h,o-z)
      dimension a(3,ncent),ngrid(ncent),grid(4,npmx,ncent),
     $        d(ncent),r(ncent,ncent),p(npmx),wtt(ncent),rad(1)
      dimension arad(ncent,ncent)
      data pi/3.14159265358/,eta/1./
c
c compute distances between atoms
c
      do 1 i=1,ncent
      do 1 j=1,ncent
       dist=sqrt((a(1,i)-a(1,j))**2
     $           +(a(2,i)-a(2,j))**2
     $           +(a(3,i)-a(3,j))**2)
 1     r(i,j)=dist         
       do 8 i=1,ncent
          do 8 j=1,ncent
             if(irad.eq.0)then
                arad(i,j)=0.
             else
                chi=rad(i)/rad(j)
                amu=(chi-1.)/(chi+1.)
                arad(i,j)=amu/(amu*amu-1.)
                if(arad(i,j).gt..5)arad(i,j)=.5
                if(arad(i,j).lt.-.5)arad(i,j)=-.5
             endif
 8     continue
c
c loop over all atoms
c
      do 2 jatom=1,ncent
c
c skip an atom if it has no weight
c
       if(wtt(jatom).gt.1.e-8)then
       npts=ngrid(jatom)
c
c loop over points associated with this atom
c
       do 3 i=1,npts
c
c compute distances between this point and all the atoms
c
        do 30 k=1,ncent
         dist=sqrt((grid(1,i,jatom)-a(1,k))**2
     $             +(grid(2,i,jatom)-a(2,k))**2
     $             +(grid(3,i,jatom)-a(3,k))**2)
 30   d(k)=dist
c
c initialize normalization
c
        sum=0
c
c compute the voronoi transformation for all atoms
c
        do 4 j=1,ncent
         pp=1
         do 5 k=1,ncent
          if(k.ne.j)then
           amu=(d(j)-d(k))/r(j,k)
           amu=amu+arad(j,k)*(1.-amu*amu)
           arg=amu
           do 6 m=1,miter
            s=1.5*arg-.5*arg**3
            arg=s
 6          continue
          pp=pp*.5*(1.-s)
          endif
 5        continue
         sum=sum+pp
c
c save the transformation for the atom of interest
c
 4      if(j.eq.jatom)p(i)=pp
c
c normalize the transformation
c
        p(i)=p(i)/sum
 3     continue
c
c compute new weights for the atom of interest      
c
       do 7 i=1,npts
 7     grid(4,i,jatom)=grid(4,i,jatom)*p(i)
      endif
 2    continue
c
c
c  do an integral and apply condensing transformation
c  example is a yukawa potential at location a(1 to 3,i))
c
      suma=0.
      do 301 ia=1,ncent
        if(wtt(ia).gt.1.e-8)then
        npt=ngrid(ia)
        do 302 ii=1,npt
          xq=grid(1,ii,ia)
          yq=grid(2,ii,ia)
          zq=grid(3,ii,ia)
          wnode=grid(4,ii,ia)
          func = 0.
          do 401 i=1,ncent
            dist =sqrt((xq-a(1,i))**2 + (yq-a(2,i))**2 + (zq-a(3,i))**2)
  401       func = func + exp(-eta*dist)/dist
        suma = suma + func*wnode
 302  continue
      endif
301   continue
      exact=0.
      do 57 i=1,ncent
      exact=exact+4.*pi/eta**2
   57 continue
      write(6,939) exact
  939 format('       exact value of yukawa integral =',d15.8)
      write(6, 936) suma
  936 format(' integral using voronoi quadrature = ',d15.8)
      return
      end
