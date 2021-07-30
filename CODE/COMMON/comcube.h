      integer mfit, nsizex, nsizey, nsizev, nboot
      logical VELmps, lfout
      real cton, imin, hmax, noise, pixscale, thresh, vdef
      real version, wmax
      common / vcube / nsizex, nsizey, nsizev, pixscale, VELmps, mfit,
     +                 version, noise, thresh, imin, wmax, hmax, nboot,
     +                 lfout, cton, vdef
c
      character*100 datfile, infile, maskfile, outroot, outpfile
      character extn( 2, 2 )*3
      common / files / infile, datfile, outpfile, outroot, maskfile,
     +                 extn
c
      integer nsize, nv
      logical vbs
      parameter ( nsize = 200 )
      real vals( 2, nsize )
      common / HermGau / nv, vbs, vals
c
      integer nkd, nkt
      real smsig
      common / smthpars / nkt, nkd, smsig
c
      real*8 pi
      parameter ( pi = 3.141592653589793238462643d0 )
