      subroutine tabcov( aparams, afitval, lfrac, m )
c Copyright (C) 2017, Jerry Sellwood and Kristine Spekkens
c
c outputs a table in csv format of all the bootstrap values of the principal
c   parameters, and light fractions, for possible analysis of co-variances
c
c   Originally written by JAS Aug 2017
c
      include 'commons.h'
c
c calling arguments
      integer m
      real*8 afitval( m, mtot ), aparams( m, md ), lfrac( 3, m )
c
c external
      integer lnblnk
c
c local arrays
      logical ldeg( md + 4 )
      real aoffs( md + 4 )
      real, allocatable :: x( :, : )
c
c local variables
      character*10 label( md + 4 )
      character comma*1, line*256
      integer i, i_bl, i_br, i_dk, iI_0, j, k, l, nI_0, nnd
      logical lc
      real val, xmin, xmax
c
c determine number of parameters for output
      nnd = nd
      if( lphot )then
        if( lI_0 )then
          nnd = nd + 1
          nI_0 = nnd
        end if
        lc = ldisk .and. ( lnax .or. lbulge )
        if( lc )then
          if( ldisk )then
            nnd = nnd + 1
            i_dk = nnd
          end if
          if( lnax )then
            nnd = nnd + 1
            i_br = nnd
          end if
          if( lbulge )then
            nnd = nnd + 1
            i_bl = nnd
          end if
        end if
      end if
      if( nnd .le. 1 )then
        print *, 'no covariances to plot'
        return
      end if
      if( nnd .gt. 19 )then
        print *, 'No .csv file created: too many parameters',
     +                                           nd, ' +', nnd - nd
        print *, 'Problem in TABCOV, contact code administrator'
        return
      end if
c set defaults
      do i = 1, nnd
        ldeg( i ) = .false.
        aoffs( i ) = 0
      end do
c set column headers, etc
      i = 0
      if( lpa )then
        i = 1
        label( i ) = 'disk pa'
        ldeg( i ) = .true.
        aoffs( i ) = -90
      end if
      if( leps )then
        i = i + 1
        label( i ) = 'disk el'
      end if
      if( lcentre )then
        label( i + 1 ) = 'xcen'
        label( i + 2 ) = 'ycen'
        i = i + 2
      end if
      if( lvels .and. lnax .and. lphib )then
        i = i + 1
        label( i ) = 'bar pa'
        ldeg( i ) = .true.
      end if
c options for photometry
      if( lphot )then
        iI_0 = 0
        if( ldisk )iI_0 = nellip
        if( lnax )then
          iI_0 = iI_0 + nbar
          if( lphib )then
            i = i + 1
            label( i ) = 'bar pa'
            ldeg( i ) = .true.
            aoffs( i ) = -90
          end if
          if( lepsb )then
            i = i + 1
            label( i ) = 'bar el'
          end if
        end if
        if( lbulge )then
          if( lbleps )then
            i = i + 1
            label( i ) = 'bulge el'
          end if
          if( lsersn )then
            i = i + 1
            label( i ) = 'Sersic n'
          end if
          if( lr_e )then
            i = i + 1
            label( i ) = 'bulge R_e'
          end if
          i = i + 1
          label( i ) = 'bulge I_0'
          iI_0 = iI_0 + 1
        end if
c light fractions
        if( lc )then
          if( ldisk )then
            i = i + 1
            label( i ) = 'f_{disk}'
          end if
          if( lnax )then
            i = i + 1
            label( i ) = 'f_{bar}'
          end if
          if( lbulge )then
            i = i + 1
            label( i ) = 'f_{bulge}'
          end if
        end if
      end if
c simple warp model
      if( lwarp )then
c inner radius of warp
        if( lrwarp )then
          i = i + 1
          label( i ) = 'R_{warp}'
        end if
c maximum eccentricity of warp
        if( lwepsm )then
          i = i + 1
          label( i ) = 'wepsm'
        end if
c max position angle of warp
        if( lwpm )then
          i = i + 1
          label( i ) = 'wphim'
          ldeg( i ) = .true.
        end if
      end if
c check that the number of parameters set agrees with the number expected
      if( i .ne. nnd )then
        print *, 'logical error in number of parameters'
        call crash( 'TABCOV' )
      end if
c allocate space for table
      allocate ( x( nnd, nunc ) )
c work over columns and rows
      do j = 1, nnd
        xmax = -10000
        xmin = -xmax
        do i = 1, nunc
          if( lI_0 .and. j .eq. nI_0 )then
            x( j, i ) = afitval( i, iI_0 )
          else if( ldisk .and. ( j .eq. i_dk ) )then
            x( j, i ) = lfrac( 1, i )
          else if( lnax .and. ( j .eq. i_br ) )then
            x( j, i ) = lfrac( 2, i )
          else if( lbulge .and. ( j .eq. i_bl ) )then
            x( j, i ) = lfrac( 3, i )
          else
            x( j, i ) = aparams( i, j )
          end if
          if( ldeg( j ) )x( j, i ) = x( j, i ) * 180. / pi + aoffs( j )
          xmax = max( xmax, x( j, i ) )
          xmin = min( xmin, x( j, i ) )
        end do
c adjust to positive values if needed
        if( ldeg( j ) )then
          val = 0
          if( xmax .lt. 0. )then
            val = 180
            do i = 1, nunc
              x( j, i ) = x( j, i ) + val
            end do
            xmin = xmin + val
            xmax = xmax + val
          end if
        end if
      end do
c create .csv file
      k = lnblnk( outroot )
      open( 9, file = outroot( 1:k )//'csv', form = 'formatted',
     +      status = 'unknown', iostat = i )
      if( i .ne. 0 )then
        print *, 'Error opening .csv file'
        call crash( 'TABCOV' )
      end if
c column header
      l = 0
      comma = ','
      do j = 1, nnd
        k = l + 1
        l = l + lnblnk( label( j ) )
        line( k:l ) = label( j )
        if( j .lt. nnd )then
          l = l + 1
          line( l:l ) = comma
        end if
      end do
      write( 9, '( a )' )line( 1:l )
c values
      do i = 1, nunc
        l = 0
        do j = 1, nnd
          k = l + 1
          l = l + 12
          write( line( k:l ), '( 1pe12.4 )' )x( j, i )
          if( j .lt. nnd )then
            l = l + 1
            line( l:l ) = comma
          end if
        end do
        write( 9, '( a )' )line( 1:l )
      end do
      close( 9 )
      return
      end
