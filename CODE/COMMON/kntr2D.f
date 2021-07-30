      subroutine kntr2D( nx, ny, direct )
c  Copyright (C) 2020, Jerry Sellwood
c
c Routine to create file(s) of the doubly transformed kernel
c
c It tabulates the kernel and calculates its Fourier transform, saving it
c   in a scratch file.
c It will also save the direct (untransformed) kernel in a separate file if
c   requested by the calling argument being set to .true.
c
c Uses single precision FFTPAK routines
      implicit none
c
      include 'comcube.h'
c
c calling argument
      integer nx, ny
      logical direct
c
c local allocatable arrays
      real, allocatable :: krn( : )
      real, allocatable :: trig( : )
      real, allocatable :: w( : )
c
c local variables
      integer i, ifail, itype, ix, iy, j, lkrn, ltrig, lwork, m
      logical pass
      real a, x, y
c
      print *, 'Creating a file of the transformed smoothing kernel'
c allocate space
      lkrn = ( nx + 1 ) * ( ny + 1 )
      allocate ( krn( lkrn ) )
      lwork = 2 * lkrn
      allocate ( w( lwork ) )
      ltrig = 3 * max( nx, ny ) + 15
      allocate ( trig( ltrig ) )
c open file(s)
      nkt = 59
      open( nkt, file = 'krnt', form = 'unformatted', status = 'new' )
      write( nkt )nx, ny, smsig
      if( direct )then
        nkd = 58
        open( nkd, file = 'krnd', form = 'unformatted', status = 'new' )
        write( nkd )nx, ny, smsig
      end if
c work over (x,y) field points
      m = 0
      do iy = 0, ny
        y = real( iy )
        do ix = 0, nx
          x = real( ix )
          m = m + 1
c compute kernel value
          a = -.5 * ( x**2 + y**2 ) / smsig**2
          krn( m ) = exp( a ) / ( 2. * pi * smsig**2 )
        end do
      end do
      if( direct )then
        write( nkd )( krn( i ), i = 1, m )
      end if
c Fourier transformation in x
      call rcosti( nx + 1, trig )
      j = 1
      do iy = 1, ny + 1
        call rcost( nx + 1, krn( j ), trig )
        j = j + nx + 1
      end do
c re-order array for Fourier transformation in y
      m = 0
      do ix = 1, nx + 1
        j = ix
        do iy = 1, ny + 1
          m = m + 1
          w( m ) = krn( j )
          j = j + nx + 1
        end do
      end do
c Fourier transformation in y
      call rcosti( ny + 1, trig )
      j = 1
      do iy = 1, nx + 1
        call rcost( ny + 1, w( j ), trig )
        j = j + ny + 1
      end do
c copy back
      m = 0
      do iy = 1, ny + 1
        j = iy
        do ix = 1, nx + 1
          m = m + 1
          krn( m ) = w( j )
          j = j + ny + 1
        end do
      end do
c write out transformed version
      m = 0
      do iy = 1, ny + 1
        j = m + 1
        m = m + nx + 1
        write( nkt )( krn( i ), i = j, m )
      end do
      deallocate ( w )
      deallocate ( trig )
      deallocate ( krn )
c close and re-open file(s) to ensure that writes are completed
      close( nkt )
      open( nkt, file = 'krnt', form = 'unformatted', status = 'old' )
      if( direct )then
        close( nkd )
        open( nkd, file = 'krnd', form = 'unformatted', status = 'old' )
      end if
      print *, 'Transformed kernel fn file ready'
      return
      end
