      subroutine smooth( dat, nx, ny, last )
c  Copyright (C) 2020, Jerry Sellwood
c
c Routine to convolve the intensities in the channel maps by a Gaussian
c   kernel.  It uses Fourier transforms on a grid that is twice the size
c   of each channel map in order to avoid wrapping smoothed signal near
c   any edge to the opposite side
c
c It computes and saves the transformed kernel to a scratch file on the
c   first call and deletes it at the end of the last call
      implicit none
c
      include 'comcube.h'
c
      integer nx, ny
      logical last
      real dat( nx, ny )
c
      real, allocatable, save :: temp( :, : )
c
      integer i, j
      real kern
      save kern
      data kern / -1 /
c
c check for this being a new smoothing kernel
      if( kern .ne. smsig )then
c compute and save the transform of the smoothing kernel
        if( smsig .le. 0. )call gtreal(
     +     'Enter sigma of Gaussian smoothing kernel in pixels', smsig )
        kern = smsig
        call kntr2D( nx, ny, .false. )
c allocate temporary space
        if( .not.
     +          allocated ( temp ) )allocate ( temp( 2 * nx, 2 * ny ) )
      end if
c copy data to work array that is 4 times larger
      do j = 1, nx
        do i = 1, ny
          temp( i, j ) = dat( i, j )
        end do
      end do
c Fourier analysis
      call anal2D( temp, nx, ny )
c convolve with transformed kernel
      call conv2D( temp, nx, ny )
c Fourier synthesis
      call synt2D( temp, nx, ny )
c copy back smoothed data
      do j = 1, nx
        do i = 1, ny
          dat( i, j ) = temp( i, j )
        end do
      end do
c
c delete the transformed kernel file and flag it as no longer available
      if( last )then
        deallocate ( temp )
        close ( nkt, status = 'delete' )
        kern = -1
      end if
      return
      end
