*deck rgauss
      function rgauss (xmean, sd)
c***begin prologue  rgauss
c***purpose  generate a normally distributed (gaussian) random number.
c***library   slatec (fnlib)
c***category  l6a14
c***type      single precision (rgauss-s)
c***keywords  fnlib, gaussian, normal, random number, special functions
c***author  fullerton, w., (lanl)
c***description
c
c generate a normally distributed random number, i.e., generate random
c numbers with a gaussian distribution.  these random numbers are not
c exceptionally good -- especially in the tails of the distribution,
c but this implementation is simple and suitable for most applications.
c see r. w. hamming, numerical methods for scientists and engineers,
c mcgraw-hill, 1962, pages 34 and 389.
c
c             input arguments --
c xmean  the mean of the guassian distribution.
c sd     the standard deviation of the guassian function
c          exp (-1/2 * (x-xmean)**2 / sd**2)
c
c***references  (none)
c***routines called  rand
c***revision history  (yymmdd)
c   770401  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   910819  added external statement for rand due to problem on ibm
c           rs 6000.  (wrb)
c***end prologue  rgauss
      external rand
c***first executable statement  rgauss
      rgauss = -6.0
      do 10 i=1,12
        rgauss = rgauss + rand(0.0)
 10   continue
c
      rgauss = xmean + sd*rgauss
c
      return
      end
