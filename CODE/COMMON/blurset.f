      subroutine blurset
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c This subroutine sets up the parameters for the seeing correction as
c   described in the documentation
c
c   Created by JAS March 2011
c   Input seeing disk size change to FWHM - Sep 2011
c   Limit range of seeing to dsee =< 3 pixels - Jun 2015
c
c common blocks
      include 'commons.h'
c
c local arrays
      integer nsprd
      parameter ( nsprd = 11 )
      integer mwide( nsprd ), npt( nsprd )
      real rmx( nsprd )
c
c local variables
      integer i, isprd, j, k
      logical firstc
      real a, b, sgma, totw
      save rmx
c
      data firstc / .true. /
      data mwide / 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9 /
      data npt / 5, 9, 13, 21, 25, 29, 37, 45, 49, 57, 61 /
      data rmx / 1., 2., 4., 5., 8., 9., 10., 13., 16., 17., 18. /
c
c initialize on first call only
      if( firstc )then
        firstc = .false.
        do i = 1, nsprd
          rmx( i ) = sqrt( rmx( i ) )
        end do
c set pixel offset arrays - useful also for bootstrap code
        k = 1
        ibl( k ) = 0
        jbl( k ) = 0
c first nearest neighbours
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
          print *, k, mblur, ' array too small'
          call crash( 'blurset' )
        end if
      end if
      if( lseeing )then
c dsee is FWHM - determine sgma for the equivalent Gaussian
        sgma = dsee / 3.33
c set number of pixels over which to spread the light
c    include pixels within twice the seeing radius
        isprd = nsprd
        do while ( ( isprd .gt. 1 ) .and.
     +             ( rmx( isprd ) .gt. dsee ) )
          isprd = isprd - 1
        end do
        if( isprd .gt. 9 )then
          print *, 'The seeing disk diameter you have entered,', dsee,
     +      'pixels, is quite large'
          print *, 'and correcting for it would be too cumbersome'
          print *, 'Your image is oversampled and could be binned'
          call crash( 'blurset' )
        end if
        iblur = mwide( isprd )
        jblur = iblur
        nblur = npt( isprd )
c set the blurring weights
        totw = 0
        do i = 1, nblur
          a = ibl( i )
          b = jbl( i )
c assume a 2D Gaussian profile
          wblur( i ) = exp( -.25 * ( a**2 + b**2 ) / sgma**2 )
          totw = totw + wblur( i )
        end do
c normalize
        do i = 1, nblur
          wblur( i ) = wblur( i ) / totw
        end do
      end if
      return
      end
