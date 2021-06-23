      subroutine revisout
c to update the output file with the uncertainties from post-processing
c of bootstrap outputc
c
c created by JAS - June 2012
c
      include 'commons.h'
c
c external
      integer lnblnk
c
c local variables
      character line*150
      integer i, j
      real a
c
      if( .not. lphot )then
        print *, 'routine programmed for photoetry data only'
        call crash( 'revisout' )
      end if
c open output parameter file
c      open( 3, file = outpfile, status = 'old', readonly )
      open( 3, file = outpfile, status = 'old' )
      read( 3, '( a)' )line
c create new file
      j = lnblnk( outroot )
      open( 4, file = outroot( 1:j )//'new', status = 'unknown' )
c copy preamble
      do while ( line( 1:19 ) .ne. 'Best fitting values' )
        j = lnblnk( line )
        write( 4, '( a)' )line( 1:j )
        read( 3, '( a)', iostat = i )line
        if( i .ne. 0 )then
          print *,
     +            'Best fiting values line not found in file ', outpfile
          call crash( 'revisout' )
        end if
      end do
c copy this line and next
      j = lnblnk( line )
      write( 4, '( a)' )line( 1:j )
      read( 3, '( a)' )line
      j = lnblnk( line )
      write( 4, '( a)' )line( 1:j )
c read and edit parameters and fitted values
      read( 3, '( a)', iostat = i )line
      do while ( i .eq. 0 )
        if( line( 1:7 ) .eq. 'disk PA' )then
          write( line( 49:53 ), '( f5.2 )' )enewpa * 180. / pi
        else if( line( 1:8 ) .eq. 'disk eps' )then
          write( line( 49:53 ), '( f5.2 )' )eneweps
        else if( line( 1:9 ) .eq. 'disk incl' )then
          read( line( 39:43 ), '( f5.2 )' )a
          a = a * pi / 180
          write( line( 49:53 ), '( f5.2 )' )
     +                                    eneweps / sin( a ) * 180. / pi
        else if( line( 1:10 ) .eq. 'x,y center' )then
          write( line( 49:53 ), '( f5.2 )' )enewxcen
          write( line( 68:72 ), '( f5.2 )' )enewycen
        else if( line( 1:11 ) .eq. 'Non-axisymm' )then
          write( line( 49:53 ), '( f5.2 )' )enewbar_pa * 180. / pi
        else if( line( 1:8 ) .eq. 'Bar eps:' )then
          write( line( 49:53 ), '( f5.2 )' )enewbar_eps
        else if( line( 1:9 ) .eq. 'Bulge r_e' )then
          write( line( 49:53 ), '( f5.2 )' )enewr_bulge
        else if( line( 1:9 ) .eq. 'Bulge I_0' )then
          write( line( 49:53 ), '( f5.2 )' )eibulge( 1 )
        else if( line( 1:10 ) .eq. 'Bulge eps:' )then
          write( line( 49:53 ), '( f5.2 )' )enewbulge_l
        else if( line( 1:19 ) .eq. 'Disk light fraction' )then
          write( line( 49:53 ), '( f5.2 )' )100. * edskfrac
        else if( line( 1:18 ) .eq. 'Bar light fraction' )then
          write( line( 49:53 ), '( f5.2 )' )100. * ebarfrac
        else if( line( 1:20 ) .eq. 'Bulge light fraction' )then
          write( line( 49:53 ), '( f5.2 )' )100. * eblgfrac
        end if
c output edited line
        j = lnblnk( line )
        if( j .gt. 0 )then
          write( 4, '( a)' )line( 1:j )
        else
          write( 4, * )
        end if
c read next line
        read( 3, '( a)', iostat = i )line
c stop when fitted intensities line is encountered
        if( line( 1:18 ) .eq. 'Fitted intensities' )i = 1
      end do
c copy this line and 2 more
      j = lnblnk( line )
      write( 4, '( a)' )line( 1:j )
      do i = 1, 2
        read( 3, '( a)' )line
        j = lnblnk( line )
        write( 4, '( a)' )line( 1:j )
      end do
c work over ellipses
      do i = 1, nellip
        read( 3, '( a)' )line
c edit in intensity uncertainties line by line
        write( line( 33:44 ), '( 6f12.2 )' )eid( i )
        write( line( 57:68 ), '( 6f12.2 )' )eibar( i )
        write( line( 81:92 ), '( 6f12.2 )' )eblgprof( i )
        j = lnblnk( line )
        write( 4, '( a)' )line( 1:j )
      end do
      close( 4 )
      return
      end
