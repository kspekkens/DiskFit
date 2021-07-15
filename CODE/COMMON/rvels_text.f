      subroutine rvels_text
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c Reads in a dataset to xval, yval, sdat, sdate from a simple text file
c
c   Created by KS
c   Polished by JAS October 2009
c   Convert data to pixel map by JAS March 2011
c   Conversion made a user option JAS July 2012
c   Updated to f90 - JAS Jan 2015
c   Fixed bug in rasterization - JAS Aug 2015
c
      include 'commons.h'
      real, allocatable :: temp( :, : )
c
c external
      integer lnblnk
c
c local variables
      character*1 yn
      integer i, j, k
      real d, dmin, spacing, xmax, xmin, ymax, ymin
c
      open( 4, file = datfile, status = 'old' )
c skip header lines
      do i = 1, 4
        read( 4, * )
      end do
c count lines of data
      i = 0
      inp_pts = -1
      do while ( i .eq. 0 )
        read( 4, *, iostat = i )d
        inp_pts = inp_pts + 1
      end do
c rewind and skip header lines
      rewind 4
      do i = 1, 4
        read( 4, * )
      end do
      xmin = 1.e6
      xmax = -xmin
      ymin = 1.e6
      ymax = -ymin
      allocate ( temp( 4, inp_pts ) )
c read in file
      do i = 1, inp_pts
        read( 4, * )( temp( j, i ), j = 1, 4 )
        xmin = min( xmin, temp( 1, i ) )
        xmax = max( xmax, temp( 1, i ) )
        ymin = min( ymin, temp( 2, i ) )
        ymax = max( ymax, temp( 2, i ) )
        if( VELmps )then
          temp( 3, i ) = 1.e-3 * temp( 3, i )
          temp( 4, i ) = 1.e-3 * temp( 4, i )
        end if
      end do
      close( 4 )
      i = lnblnk( datfile )
      print *, inp_pts,
     +                ' velocity points read from file ', datfile( 1:i )
c
c test whether conversion to a rectangular pixel map is possible
c
c find smallest non-zero x & y differences in input values
      spacing = 1.e6
      dmin = .0005 * min( xmax - xmin, ymax - ymin )
      do i = 1, inp_pts - 1
        do j = i + 1, inp_pts
          d = abs( temp( 1, i ) - temp( 1, j )  )
          if( d .gt. dmin )spacing = min( spacing, d )
          d = abs( temp( 2, i ) - temp( 2, j )  )
          if( d .gt. dmin )spacing = min( spacing, d )
        end do
      end do
c check image size
      lox = 1
      loy = 1
      xrange = int( ( xmax - xmin ) / spacing ) + 1
      yrange = int( ( ymax - ymin ) / spacing ) + 1
      l2D = .false.
c value of 1024 here is arbitrary
      if( ( xrange .le. 1024 ) .and. ( yrange .le. 1024 ) )then
c 2D raster option possible if the user wishes
        yn = 'x'
        do while ( ( yn .ne. 'y' ) .and. ( yn .ne. 'n' ) )
          print *, 'Your input data points are spaced regularly' //
     + ' enough that they could be inserted into a 2D raster of pixels.'
          print *, 'This has both advantages and disadvantages - see' //
     +        ' section 2.1.2 of the documentation'
          print *, 'Would you like DiskFit to rasterize your data' //
     +' (y/n)?  (Enter n if you are inexperienced or do not understand)'
          read '( a )', yn
c ensure lowercase
          j = ichar( yn )
          if( ( j .gt. 64 ) .and. ( j .lt. 91 ) )then
            yn = char( j + 32 )
          end if
        end do
        l2D = yn .eq. 'y'
      end if
      if( l2D )then
c revise xrange, yrange to ensure all pixels are covered
        j = nint( xmin * spacing )
        xmin = real( j ) / spacing
        xrange = int( ( xmax - xmin ) / spacing ) + 2
        j = nint( ymin * spacing )
        ymin = real( j ) / spacing
        yrange = int( ( ymax - ymin ) / spacing ) + 2
c dimension arrays
        allocate ( xval( xrange ) )
        allocate ( yval( yrange ) )
        allocate ( sdat( xrange, yrange ) )
        allocate ( sdate( xrange, yrange ) )
        allocate ( lgpix( xrange, yrange ) )
c set (x,y) positions of pixels
        do i = 1, xrange
          xval( i ) = real( i - 1 ) * spacing + xmin
        end do
        do i = 1, yrange
          yval( i ) = real( i - 1 ) * spacing + ymin
        end do
c preset arrays
        do j = 1, yrange
          do i = 1, xrange
            sdat( i, j ) = -1.e6
            lgpix( i, j ) = .false.
          end do
        end do
c insert input values into pixel map
        do k = 1, inp_pts
          i = nint( ( temp( 1, k ) - xmin ) / spacing ) + 1
          j = nint( ( temp( 2, k ) - ymin ) / spacing ) + 1
          lgpix( i, j ) = .true.
          sdat( i, j ) = temp( 3, k )
          sdate( i, j ) = temp( 4, k )
        end do
c parameter needed for plotting
        istepout = max( nint( spacing ), 1 )
      else
c allocate space as a 1D array
        allocate ( xval( inp_pts ) )
        allocate ( yval( inp_pts ) )
        allocate ( sdat( inp_pts, 1 ) )
        allocate ( sdate( inp_pts, 1 ) )
        do i = 1, inp_pts
          xval( i ) = temp( 1, i )
          yval( i ) = temp( 2, i )
          sdat( i, 1 ) = temp( 3, i )
          sdate( i, 1 ) = temp( 4, i )
        end do
c force program crash if not l2D and either seeing
c     correction, or junc<0 for bootstrap selected
        if( lseeing )then
          print*, 'No seeing correction allowed for pixel lists'
          print*, 'Set seeing FWHM = 0 and re-run'
          call crash( 'rvels_text' )
        end if
        if( luncert .and. ( junc .lt. 0 ) )then
          print*, 'Chosen bootstrap not allowed for pixel lists'
          print*, 'Set junc > 0 and re-run'
          call crash( 'rvels_text' )
        end if
      end if
      deallocate ( temp )
      return
      end
