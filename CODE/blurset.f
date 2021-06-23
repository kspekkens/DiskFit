      subroutine blurset
c This subroutine sets up the parameters for the seeing correction as
c   described in the documentation
c
c   Created by JAS March 2011
c   Input seeing disk size change to FWHM - Sep 2011
c
c common blocks
      include 'commons.h'
c
c local variables
      integer i, j, k
      real a, b, totw
c
c set pixel offset arrays - useful also for bootstrap code
      do i = 1, mblur
        ibl( i ) = 0
        jbl( i ) = 0
      end do
c first nearest neighbours
      k = 1
      do i = -1, 1, 2
        ibl( k + 1 ) = i
        jbl( k + 2 ) = i
        k = k + 2
      end do
c second nearest neighbours
      do i = -1, 1, 2
        do j = -1, 1, 2
          k = k + 1
          ibl( k ) = i
          jbl( k ) = j
        end do
      end do
c third nearest neighbours
      do i = -1, 1, 2
        ibl( k + 1 ) = 2 * i
        jbl( k + 2 ) = 2 * i
        k = k + 2
      end do
c fourth nearest neighbours
      do i = -1, 1, 2
        do j = -1, 1, 2
          ibl( k + 1 ) = 2 * i
          jbl( k + 1 ) = j
          ibl( k + 2 ) = i
          jbl( k + 2 ) = 2 * j
          k = k + 2
        end do
      end do
c fifth nearest neighbours
      do i = -1, 1, 2
        do j = -1, 1, 2
          k = k + 1
          ibl( k ) = 2 * i
          jbl( k ) = 2 * j
        end do
      end do
c sixth nearest neighbours
      do i = -1, 1, 2
        ibl( k + 1 ) = 3 * i
        jbl( k + 2 ) = 3 * i
        k = k + 2
      end do
c seventh nearest neighbours
      do i = -1, 1, 2
        do j = -1, 1, 2
          ibl( k + 1 ) = 3 * i
          jbl( k + 1 ) = j
          ibl( k + 2 ) = i
          jbl( k + 2 ) = 3 * j
          k = k + 2
        end do
      end do
c eighth nearest neighbours
      do i = -1, 1, 2
        do j = -1, 1, 2
          ibl( k + 1 ) = 3 * i
          jbl( k + 1 ) = 2 * j
          ibl( k + 2 ) = 2 * i
          jbl( k + 2 ) = 3 * j
          k = k + 2
        end do
      end do
c ninth nearest neighbours
      do i = -1, 1, 2
        ibl( k + 1 ) = 4 * i
        jbl( k + 2 ) = 4 * i
        k = k + 2
      end do
c tenth nearest neighbours
      do i = -1, 1, 2
        do j = -1, 1, 2
          ibl( k + 1 ) = 4 * i
          jbl( k + 1 ) = j
          ibl( k + 2 ) = i
          jbl( k + 2 ) = 4 * j
          k = k + 2
        end do
      end do
c eleventh nearest neighbours
      do i = -1, 1, 2
        do j = -1, 1, 2
          k = k + 1
          ibl( k ) = 3 * i
          jbl( k ) = 3 * j
        end do
      end do
      if( k .gt. mblur )then
        print *, mblur, ' array too small'
        call crash( 'blurset' )
      end if
c rsee now FWHM - decide whether seeing corrections are worthwhile
      if( rsee .le. 0.67 )then
        lseeing = .false.
      else
c set number of pixels over which to spread the light
        if( rsee .gt. 2.8 )then
          iblur = 5
          nblur = 25
        else if( rsee .gt. 2.4 )then
          iblur = 5
          nblur = 21
        else if( rsee .gt. 1.33 )then
          iblur = 5
          nblur = 13
        else if( rsee .gt. .94 )then
          iblur = 3
          nblur = 9
        else if( rsee .ge. .67 )then
          iblur = 3
          nblur = 5
        end if
        jblur = iblur
c set the blurring weights
        totw = 0
        do i = 1, nblur
          a = ibl( i )
          b = jbl( i )
          wblur( i ) = exp( -.25 * ( a**2 + b**2 ) / rsee**2 )
          totw = totw + wblur( i )
        end do
c normalize
        do i = 1, nblur
          wblur( i ) = wblur( i ) / totw
        end do
      end if
      return
      end
