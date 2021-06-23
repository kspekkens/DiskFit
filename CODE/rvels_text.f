      subroutine rvels_text
c Reads in a dataset to xval, yval, sdat, sdate from a simple text file
c
c   Created by KS
c   Polished by JAS October 2009
c   Convert data to pixel map by JAS March 2011
c   Conversion made a user option JAS July 2012
      include 'commons.h'
      real temp( 4, mpx )
      equivalence ( temp( 1, 1 ), ldat( 1, 1 ) )

c local variables
      character*1 yn
      integer i, j, k, lnblnk
      real d, spacing, xmax, xmin, ymax, ymin
c
      open( 4, file = invfile, status = 'old' )
c skip header lines
      do i = 1, 4
        read( 4, * )
      end do
      xmin = 1.e6
      xmax = -xmin
      ymin = 1.e6
      ymax = -ymin
c read in file
      do i = 1, mpx * mpy / 4
        read( 4, *, end = 99 )( temp( j, i ), j = 1, 4 )
        xmin = min( xmin, temp( 1, i ) )
        xmax = max( xmax, temp( 1, i ) )
        ymin = min( ymin, temp( 2, i ) )
        ymax = max( ymax, temp( 2, i ) )
        if( VELmps )then
          temp( 3, i ) = 1.e-3 * temp( 3, i )
          temp( 4, i ) = 1.e-3 * temp( 4, i )
        end if
      end do
c no end of file encountered
      print *
      print *, '( ''WARNING: maximum number '', i8, ' //
     +           ' ''of input velocities read'' )', mapx, mapy
      j = lnblnk( invfile )
      print *, '( ''It is possible that some lines of '', a, ' //
     +                     ' '' have not been read.'' )', invfile( 1:j )
      print *, 'Check outputs carefully!'
      print *
c normal end
   99 close( 4 )
      inp_pts = i - 1
      i = lnblnk( invfile )
      print *, inp_pts,
     +                ' velocity points read from file ', invfile( 1:i )
c
c convert to a rectangular pixel map
c
c find smallest non-zero x & y differences in input values
      spacing = 1.e6
      do i = 1, inp_pts - 1
        do j = i + 1, inp_pts
          d = abs( temp( 1, i ) - temp( 1, j )  )
          if( d .gt. 0.1 )spacing = min( spacing, d )
          d = abs( temp( 2, i ) - temp( 2, j )  )
          if( d .gt. 0.1 )spacing = min( spacing, d )
        end do
      end do
c check image size
      lox = 1
      loy = 1
      xrange = int( ( xmax - xmin ) / spacing ) + 1
      yrange = int( ( ymax - ymin ) / spacing ) + 1
      l2D = .false.
c 2D raster option possible if the user wishes
      if( ( xrange .le. mapx ) .and. ( yrange .le. mapy ) )then
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
c set (x,y) positions of pixels
      if( l2D )then
        j = nint( xmin * spacing )
        xmin = real( j ) / spacing
        do i = 1, xrange
          xval( i ) = real( i - 1 ) * spacing + xmin
        end do
        j = nint( ymin * spacing )
        ymin = real( j ) / spacing
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
c check space
        if( inp_pts .gt. mval * min( mpx, mpy ) )then
          print *, 'xval and yval arrays too small - increase mval'
          print *, 'contact code administrator'
          call crash( 'rvels_text' )
        end if
        do i = 1, inp_pts
          xval( i ) = temp( 1, i )
          yval( i ) = temp( 2, i )
          sdat1( i ) = temp( 3, i )
          sdate1( i ) = temp( 4, i )
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
      return
      end
