      subroutine screenoutput
c screen output with input to minimization
c
c   Orignally written by KS
c   Adapted by JAS March 2011 for photometric option
c   Revised by JAS July 2011 to include warp analysis
c
      include 'commons.h'
c
c local variables
      integer i, s1l, s2l, lnblnk
      character*80 str
c
      write( *, * )
      write( *, * )'Read in input parameters from: ', infile( 1:100 )
      write( *, * )
      if( ldisk )then
        write( *, '( a22, f5.2 )' )'Input disk eps: (deg) ', eps
        write( *, '( a28, f6.2 )' )'Input disk PA (N->E; deg) = ',
     +                                              180. * pa / pi - 90.
      end if
      write( *, '( a30, f8.2 )' )
     +                     'Input x-center (data units) = ', xcen
      write( *, '( a30, f8.2 )' )
     +                     'Input y-center (data units) = ', ycen
      if( lvels )write( *, '( a36, f8.2 )' )
     +              'Input disk systemic velocity (km/s) = ', vsys
      if( lphot )write( *, '( a27, 3f8.2 )' )
     +              'Sky, sky sig, gain (ADU) = ', sky, skysig, gain
      write( *, * )
      if( lpa .or. leps .or. lcentre .or. lnax .or. lbulge .or.
     +    lsystemic .or. lwarp )then
        write( *, * )'Fitting for:'
        if( lpa )write( *, * )'disk PA'
        if( leps )write( *, * )'disk ellipticity'
        if( lcentre )write( *, * )'disk x-center, y-center'
        if( lnax )then
           if( lvels )write( *, * )'non-axisymmetric disk components'
           if( lphot )write( *, * )'Bar profile'
        end if
        if( lphib )then
           if( lvels )write( *, * )'PA of non-axisymmetric component'
           if( lphot )write( *, * )'Bar PA'
        end if
        if( lepsb )write( *, * )'Bar ellipticity'
        if( lradial )write( *,* )'radial flows'
        if( lbulge )write( *, * )'Sersic bulge profile'
        if( lr_e )write( *, * )'bulge effective radius'
        if( lbleps )write( *, * )'bulge ellipticity'
        if( lsersn )write( *, * )'bulge Sersic index'
        if( lsystemic )write( *, * )'disk systemic velocity'
        if( lwarp )write( *, * )'warp model parameters'
        write( *, * )'--> Other attributes fixed to input values'
      else
        if( lvels )write( *, * )
     +    'Disk PA, inc, center & systemic vely held fixed with no warp'
        if( lphot )write( *, * )
     +       'Disk PA, inc, center, bar and bulge parameters held fixed'
      end if
      write( *, * )
      if( lvels )write( *, * )
     +                   'Velocity field components extracted at radii:'
      if( lphot )write( *, * )
     +                        'Surface brightnesses extracted at radii:'
      write( *, '( 10f7.1 )' )( sma( i ), i = 1, nellip )
      if( lradial )then
        write( *, * )
        write( *, * )
     +  'Fitting for m=0 radial flow in addition to rotation velocities'
        if( minr .gt. sma( 1 ) )write( *, * )
     +    'Radial flows set to 0 in rings inside ', sma( nminr )
        if( maxr .lt. sma( nellip ) )write( *, * )
     +    'Radial flows set to 0 in rings beyond ', sma( nmaxr )
      end if
      if( lnax )then
        write( *, * )
        if( lvels ) then
         if( order .eq. 2 )then
           write( *, * )'Fitting for m = 2 bisymmetric flow' //
     +                             ' (radial and tangential components)'
         else
           write( *, '( ''Fitting for m ='', i2, '' non-axisymmetric'' '
     +      //  ' '' flow (radial and tangential components)'' )' )order
         end if
         write( *, * )'in addition to rotation velocities'
c         write( *, * )
         if( minr .gt. sma( 1 ) )write( *, * )
     +       'Non-axisymm flows set to 0 in rings inside ', sma( nminr )
         if( maxr .lt. sma( nellip ) )write( *, * )
     +       'Non-axisymm flows set to 0 in rings beyond ', sma( nmaxr )
        else
c         write( *, * )
         write( *, * )'Fitting for a bar component'
         if( minr .gt. sma( 1 ) )write( *, * )
     +               'Axisymmetric model  inside ', sma( nminr )
         if( maxr .lt. sma( nellip ) )write( *, * )
     +               'Axisymmetric model outside ', sma( nmaxr )
         end if
      end if
c      if( lbulge )then
c        write( *, * )
c construct output string
c        write( str, '( a )' )'Bulge fitting for:'
c        if( lbleps )then
c          i = lnblnk( str )
c          write( str( i+2:i+12 ), '( a )' )'ellipticity'
c        end if
c        if( lsersn )then
c          i = lnblnk( str )
c          write( str( i+2:i+13 ), '( a )' )'Sersic index'
c        end if
c        if( lr_e )then
c          i = lnblnk( str )
c          write( str( i+2:i+17 ), '( a )' )'effective radius'
c        end if
c        i = lnblnk( str )
c        write( *, * )str( 1:i )
c      end if
      if( lwarp )then
        write( *, * )
c construct output string
        write( str, '( a )' )'Warp fitting for:'
        if( lrwarp )then
          i = lnblnk( str )
          write( str( i+2:i+13 ), '( a )' )'inner radius'
        end if
        if( lwepsm )then
          i = lnblnk( str )
          write( str( i+2:i+8 ), '( a )' )'max eps'
        end if
        if( lwepsm )then
          i = lnblnk( str )
          write( str( i+2:i+7 ), '( a )' )'max PA'
        end if
        i = lnblnk( str )
        write( *, * )str( 1:i )
      end if
      write( *, * )
      if( lvels )then
        write( *, * )
     +             'Read in measured velocities and uncertainties from:'
        s1l = lnblnk( invfile )
        s2l = lnblnk( inevfile )
        write( *, * )invfile( 1:s1l )
        if( evelfile )then
          write( *, * )inevfile( 1:s2l )
        else
          write( *, * )
        end if
        if ( .not. l2D )then
         write( *, * ) 'Input velocity list not convenient for'
         write( *, * ) '2D representation: pixel list used'
        end if
      end if
      if( lphot )then
        write( *, * )'Read in photometric image from:'
        s1l = lnblnk( inpfile )
        write( *, * )inpfile( 1:s1l )
      end if
      write( *, *)
      if( lvels )then
        write( *, * )'Number of datapoints read in:', inp_pts
        write( *, * )'Number of datapoints to be used in fit:', pixct
      else
        write( *, * )
     +        'Number of pixels in selected part of the image:', inp_pts
        write( *, * )'Number of pixels used in fit:', pixct
      end if
      return
      end
