      subroutine conv2D( temp, nx, ny )
c  Copyright (C) 2020, Jerry Sellwood
c
c Multiplies the transformed data array by the appropriate kernel
c   coefficients
c The transformed kernel is presumed stored in a file (created by KNTR2D)
c   and is read in to be used here.
c The results over-write the input array
      implicit none
c
      include 'comcube.h'
c
c calling arguments
      integer nx, ny
      real temp( 4 * nx * ny )
c
c local allocatable array
      real, allocatable, save :: trk(:)
c
c local variables
      integer i, ifail, ix, iy, iuse, j, k
      logical nyskip
      real a
      save iuse
c
      data iuse / 0 /
c
c check Green function file header
      if( iuse .eq. 0 )then
        read( nkt )i, j, a
        if( ( i .ne. nx ) .or. ( j .ne. ny ) .or.
     +      ( a .ne. smsig ) )call crash( 'CONV2D',
     +                                 'Wrong transformed kernel file' )
        print *, 'Kernel table used from file'
        iuse = 1
        allocate ( trk( nx + 1 ) )
      end if
c restart and skip header
      rewind nkt
      read( nkt )
c work over y Fourier components
      do iy = 1, ny + 1
c get Green function
        read( nkt )( trk( ix ), ix = 1, nx + 1 )
c i & k point to the cc & cs terms respectively, i + 1 & k + 1 to sc & ss
        i = max( 1, 2 * ( iy - 1 ) )
        i = ( i - 1 ) * 2 * nx + 1
        k = i + 2 * nx
c work over x Fourier components
        nyskip = ( iy .gt. 1 ) .and. ( iy .le. ny )
        do ix = 1, nx + 1
          temp( i ) = temp( i ) * trk( ix )
          if( nyskip )temp( k ) = temp( k ) * trk( ix )
          if( ( ix .gt. 1 ) .and. ( ix .le. nx ) )then
            temp( i + 1 ) = temp( i + 1 ) * trk( ix )
            if( nyskip )temp( k + 1 ) = temp( k + 1 ) * trk( ix )
            i = i + 2
            k = k + 2
          else
            i = i + 1
            k = k + 1
          end if
        end do
      end do
      return
      end
