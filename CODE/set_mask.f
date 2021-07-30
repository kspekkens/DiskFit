      program set_mask
c Copyright (C) 2020, Jerry Sellwood
c
c Main program for making a mask map from a THINGS data cube
c
c This software package was developed Jerry Sellwood in collaboration with
c    Kristine Spekkens
c
c Released May 2020
      use aacube
      implicit none
c
c common block
      include 'comcube.h'
c
c externals
      integer lnblnk
      real mask
c
c local arrays
      real, allocatable :: work( :, : ), work2( :, : )
c
c local variables
      character*80 fname
      integer i, j, k, l
      real amean, clip, fwhm, sigm, x
c
      call gtchar( 'Enter name of .inp file', fname )
      i = lnblnk( fname )
      if( fname( i-2:i ) .ne. 'inp' )then
        fname( i+1:i+4 ) = '.inp'
        i = i + 4
      end if
      call getinp( fname( 1:i ) )
c
c read the data cube
      call rcube_FITS
      print *, 'read data cube', nsizex, nsizey, nsizev
c smooth the channels one at a time
c      fwhm = 30
      call gtreal( 'Enter FWHM for smoothed data in arcsec', fwhm )
c convert to pixels
      fwhm = fwhm / pixscale
c convert to Gaussian sigma
      smsig = fwhm / 2.3548200450309489
      print *, 'Gaussian smoothing width in pixels =', smsig
c smooth each channel in turn
      do i = 1, nsizev
        call smooth( dcube( 1, 1, i ), nsizex, nsizey, i .eq. nsizev )
        if( mod( i, 10 ) .eq. 1 )print *, 'done channel', i
      end do
c copy 2 planes from all 6 edges of the cube from which to estimate the noise
      l = nsizex * nsizey * 4
      l = l + 4 * nsizex * ( nsizev - 4 )
      l = l + 4 * ( nsizey - 4 ) * ( nsizev - 4 )
      allocate( work( l, 1 ) )
      call nsevls( work, l )
c use the bi-weight estimator to set the noise value
      call biwght( work, l, amean, sigm )
      print *, 'mean and sigma from biwght', amean, sigm
      deallocate ( work )
c
      noise = sigm
c set threshold
c      imin = 2
      call gtreal( 'Enter min I/N for 3 consecutive channels', Imin )
      allocate ( work( nsizex, nsizey ) )
      allocate ( work2( nsizex, nsizey ) )
c create the initial mask, as described in Walter+08, with a copy in work2
      k = 0
      do j = 1, nsizey
        do i = 1, nsizex
          work( i, j ) = mask( i, j )
          if( work( i, j ) .gt. 0. )then
            k = k + 1
          end if
          work2( i, j ) = work( i, j )
        end do
      end do
      print *, k, 'pixels within initial mask'
c finished with dcube
      deallocate ( dcube )
c aggressively smooth the initial mask
c      smsig = 30
      call gtreal( 'Enter smoothing for super mask, in pixels', smsig )
      call smooth( work2, nsizex, nsizey, .true. )
c create the super mask and use it to clip the initial mask
      k = 0
c      clip = 3
      call gtreal( 'Enter clipping value for super mask', clip )
      do j = 1, nsizey
        do i = 1, nsizex
          if( work2( i, j ) .gt. clip )then
            if( work( i, j ) .gt. 0 )k = k + 1
          else
            work( i, j ) = 0
            work2( i, j ) = 0
          end if
        end do
      end do
      print *, k, 'pixels within final mask'
      deallocate ( work2 )
      call write_mask( work, nsizex, nsizey )
      end
