      program makemap
c Copyright (C) 2020, Jerry Sellwood
c
c Main program for making velocity maps from data cubes
c
c This software package was developed Jerry Sellwood in collaboration with
c    Kristine Spekkens
c
c Released Jan 2020
      use aacube
      implicit none
c
c common block
      include 'comcube.h'
c
c external
      integer lnblnk
c
c local arrays
      character*6 string( 6 )
      real pars( 8, 2 )
      real, allocatable :: work( :, :, : )
c
c local variables
      character*80 fname
      integer i, iOK, iwork, j, jd, k, l, m, ifail( 7 )
      logical lfit, lmaskl
      real amean, sigm, vfac, x, y

      real fmax, fmin
      integer ix, iy, kk
      common / pix / ix, iy, kk

c
      data ifail / 7 * 0 /
      data string / 'cont/N', 'S/N   ', 'mean V',
     +              'sigma ', 'h3    ', 'h4    ' /
c
c      fname = 'N2841'
      call gtchar( 'Enter name of input file', fname )
      i = lnblnk( fname )
      if( fname( i-2:i ) .ne. 'inp' )then
        fname( i+1:i+4 ) = '.inp'
        i = i + 4
      end if
      call getinp( fname( 1:i ) )
      vdef = -9999
      if( velMPS )vdef = 1.e3 * vdef
c
c read the data cube
      call rcube_FITS
      print *, 'read data cube'
      if( nboot .ge. nsizev )then
        print *, 'Number of bootstraps requested exceeds the number of'
     +          // ' velocity channels'
        call crash( 'nboot too large in MAKEMAP' )
      end if
c
c copy 2 planes from all 6 edges of the cube from which to estimate the noise
      l = nsizex * nsizey * 4
      l = l + 4 * nsizex * ( nsizev - 4 )
      l = l + 4 * ( nsizey - 4 ) * ( nsizev - 4 )
      allocate( work( l, 1, 1 ) )
      call nsevls( work, l )
c use the bi-weight estimator to set the noise value
      call biwght( work, l, amean, sigm )
      if( sigm .gt. 0. )then
        cton = amean / sigm
        print *, 'mean and sigma from biwght', cton, sigm
      else
        sigm = 3.6e-4
        cton = 0
        print *, 'Default sigma set to', sigm
      end if
      deallocate ( work )
c
      noise = sigm
c rescale vels
      vfac = 1.e-2
      if( VELmps )vfac = 1.e-5
      do i = 1, nsizev
        vval( i ) = vfac * vval( i )
      end do
      vfac = 1. / vfac
c allocate space for maps
      if( maskfile( 1:1 ) .ne. ' ' )then
        lmaskl = .true.
        allocate ( llmask( nsizex, nsizey ) )
        call get_mask
      else
        lmaskl = .false.
      end if
      if( nboot .gt. 0 )then
        allocate ( work( nsizex, nsizey, 2 * mfit + 3 ) )
      else
        allocate ( work( nsizex, nsizey, mfit + 2 ) )
      end if
      print *, 'Making map'
c work over map
      y = 1.e10
      x = -y
      k = 0
      l = 0
      jd = 0
      vbs = .false.
      do iy = 1, nsizey
        do ix = 1, nsizex
c set default values
          do m = 1, mfit + 2
            pars( m, 1 ) = 0
            pars( m, 2 ) = 0
          end do
          pars( mfit + 1, 1 ) = vdef
c check mask
          if( lmaskl )then
            lfit = llmask( ix, iy )
          else
            lfit = .true.
          end if
          if( lfit )then
            l = l + 1
            if( vbs )print *, ix, iy
            call vGHfit( ix, iy, pars, ifail )
          end if
c save parameters fitted by sumsl, Vp, and chisq, and uncertainties if created
          do m = 1, mfit + 2
            work( ix, iy, m ) = pars( m, 1 )
            if( ( nboot .gt. 0 ) .and.
     +  ( m .lt. mfit + 2 ) )work( ix, iy, m + mfit + 2 ) = pars( m, 2 )
          end do
          if( pars( mfit + 1, 1 ) .gt. -9998. )then
            k = k + 1
c convert to S/N
            work( ix, iy, 2 ) = pars( 2, 1 ) /
     +                                ( sqrt( 2. * pi ) * pars( 4, 1 ) )
c take out velocity rescaling
            work( ix, iy, 3 ) = pars( 3, 1 ) * vfac
            work( ix, iy, 4 ) = pars( 4, 1 ) * vfac
            work( ix, iy, mfit + 1 ) = pars( mfit + 1, 1 ) * vfac
c find range of fitted velocities
            x = max( x, work( ix, iy, mfit + 1 ) )
            y = min( y, work( ix, iy, mfit + 1 ) )
            if( nboot .gt. 0 )then
              work( ix, iy, mfit + 4 ) = pars( 2, 2 ) /
     +                                ( sqrt( 2. * pi ) * pars( 4, 1 ) )
              work( ix, iy, mfit + 5 ) = pars( 3, 2 ) * vfac
              work( ix, iy, mfit + 6 ) = pars( 4, 2 ) * vfac
              work( ix, iy, 2 * mfit + 3 ) = pars( mfit + 1, 2 ) * vfac
            end if
c count bi-modal line profiles
            if( mfit .gt. 4 )then
              if( nboot .eq. 0 )then
                if( pars( 2, 1 ) .lt. 0 )jd = jd + 1
              else
                if( pars( mfit + 1, 2 ) .lt. 0 )jd = jd + 1
              end if
            end if
          end if
        end do
        if( mod( iy, 100 ) .eq. 0 )print *,
     +                               'done', iy, ' rows', l, 'pixels'
      end do
      vfac = 1
      if( VELmps )vfac = 1.e-3
      print *, k, 'out of', l, ' estimated velocities in range',
     +                        vfac * y, vfac * x, ' km/s'
      if( mfit .gt. 4 )print *, jd, 'lines flagged as bi-modal'
      if( k .lt. 100 )call crash( 'MAKEMAP: nearly all zeros' )
      print '( i8, a )', ifail( 1 ), ' no signal'
      print '( i8, a )', ifail( 2 ), ' line fitter failed'
      print '( i8, a, f3.1 )', ifail( 3 ), ' S/N < ', thresh
      print '( i8, a )', ifail( 4 ), ' fitted V outside channel range'
      x = 100. * wmax * abs( vval( 1 ) - vval( 2 ) )
      print '( i8, a, f6.1, a )',
     +                    ifail( 5 ), ' line width < 0 or >', x, ' km/s'
      print '( i8, a, f5.2 )',
     +                     ifail( 6 ), ' abs values of h3 or h4 >', hmax
      if( nboot .gt. 0 )print '( i8, a, f5.2 )',
     +                      ifail( 7 ), ' nonsense values in bootstraps'
c rescale vels
      do i = 1, nsizev
        vval( i ) = vval( i ) * vfac
      end do
c output maps as FITS files
      i = mfit + 1
      if( nboot .gt. 0 )i = 2 * i
      i = i + 1
      call writemap( work, nsizex, nsizey, i )
      end
