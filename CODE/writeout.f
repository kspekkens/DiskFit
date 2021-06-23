      subroutine writeout
c Write out optimal parameters and velocity field if ltext
c Details of minimization performed at end of disk parameter file.
c
c   Originally written by KS
c   Adapted by JAS March 2011 for photometric option
c   Revised by JAS July 2011 for warp analysis
c   Revised by KS  Oct 2011 to integrate vels, phot options
c   Polished by KS Mar 2012
c
      include 'commons.h'
c
c external
      integer lnblnk
c
c local variables
      character comma*1, dash*1, fittype*4, pm*5
      integer i, j
      real a, astpa, bmodpa, dmodpa, naxphi, newi, phibsky1, phibsky2
      real planeastpa, planepa, ratio
c
      dash = '-'
      pm = ' +/- '
      comma = ','
c
c create "output model file" only if maps are in text format
      if( ltext )then
        i = lnblnk( outroot )
        outmfile = outroot( 1:i )//'mod'
        open( 3, file = outmfile, status = 'unknown' )
        write( 3, * )'     X           Y          VEL      ' //
     +              '     EVEL        SIGn        VMOD       VEL-VMOD  '
     +          // '(VEL-VMOD)/SIGn' //
     +          '   Vt           Vr          Vm,t          Vm,r        '
        write( 3, '( 149a )' )( dash, i = 1, 149 )
        a = 1
        if( VELmps )a = 1000
        if( l2D )then
          do j = 1, yrange
            do i = 1, xrange
              if( lgpix( i, j ) )then
                write( 3, '( 3f12.2, 2f11.2, 7f13.2 )' )
     +            xval( i ), yval( j ), a * sdat( i, j ),
     +            a * sdate( i, j ), a * sigma( i, j ),
     +            a * model( i, j ), a * res( i, j ),
     +            x2res( i, j ), a * rotvels( i, j ),
     +            a * radvels( i, j ), a * bivtan( i, j ),
     +            a * bivrad( i, j )
              end if
            end do
          end do
        else
          do i = 1, inp_pts
            if( lgpix1( i ) )then
              write( 3, '( 3f12.2, 2f11.2, 7f13.2 )' )
     +          xval( i ), yval( i ), a * sdat1( i ),
     +          a * sdate1( i ), a * sigma1( i ), a * model1( i ),
     +          a * res1( i ), x2res1( i ), a * rotvels1( i ),
     +          a * radvels1( i ), a * bivtan1( i ), a * bivrad1( i )
            end if
          end do
        end if
        close( 3 )
      end if
c create output parameter file
      open( 4, file = outpfile, status = 'unknown' )
      if( lvels )then
        fittype = 'vels'
      else
        fittype = 'phot'
      end if
      write( 4, '( a, a4 )' )'Minimization output, ', fittype
      write( 4, '( a, a4 )' )'-------------------------'
c write out input files
      write( 4, '( a )' )'Input files: '
      write( 4, '( a )' )infile( 1:79 )
      if( lvels )then
        if( lFITS )then
          write( 4, '( a, f5.2, a )' )invfile( 1:53 )//' pixscale ',
     +                                           pixscale, ' arcsec/pix'
          if( evelfile )then
            write( 4, '( a )' )inevfile( 1:53 )
          else
            write( 4, * )
          end if
        else
          write( 4, '( a )' )invfile( 1:53 )
          write( 4, * )
        end if
      end if
      if( lphot )then
        j = lnblnk( inpfile )
        write( 4, '( a, f5.2, a )' )inpfile( 1:j )//'  pixscale ',
     +                                           pixscale, ' arcsec/pix'
        write( 4, '( a4, f8.2, 2( a9, f8.2 ) )')
     +        'sky:',  sky, '  skysig:', skysig, '    gain:', gain
      end if
      write( 4, * )
c FITS region and sampling parameters
      if( lFITS )then
        write( 4, '( a, 4( a9, i7 ) )' )'FITS region:    ',
     +              '   x-low:',     lox,  '   y-low:', loy,
     +              ' x-range:',  xrange,  ' y-range:', yrange
c put regpa back into deg:
        regpa = ( regpa * 180. / pi ) - 90.
        write( 4, '( a, 3(a9, f7.2), a9, i7 )' )'FITS sampling:  ',
     +              '  regrad:',  regrad,  '   regpa:', regpa,
     +              '  regeps:',  regeps,  ' istpout:', istepout
      else
        write( 4, * )
        write( 4, * )
      end if
      write( 4, * )
c output files
      write( 4, '( a )' )
     +                  'Output model and (data-model) residuals files:'
      j = lnblnk( outroot )
      write( 4, '( a )' )outroot( 1:j )//'mod'
      if( lFITS )then
        write( 4, '( a )' )outroot( 1:j )//'res'
      else
        write( 4, * )
      end if
      write( 4, * )
c toggles:
      write( 4, '( a, 5( a9, l2 ) )' )'Disk toggles:   ',
     +               '      PA:', lpa,     '     eps:', leps,
     +               '  center:', lcentre, ' non-axi:', lnax,
     +               '    phib:', lphib
      if( lvels )then
        write( 4, '( a, 7( a9, l2 ) )' )'vels toggles:   ',
     +               ' interp0:', linter0, '  radial:', lradial,
     +               '    Vsys:', lsystemic,'   warp:', lwarp,
     +               '      rw:', lrwarp,  '    welm:', lwepsm,
     +               '   wphim:', lwpm
      end if
      if( lphot )then
        write( 4, '( a, 5( a9, l2 ) )' )'phot toggles:   ',
     +               '  lbulge:', lbulge,  '    lr_e:', lr_e,
     +               '  lsersn:', lsersn,  ' blg eps:', lbleps,
     +               ' bar eps:', lepsb
      end if
      write( 4, *  )
c input values
      write( 4, '( a )' )'Input values '
      write( 4, '( a )' )'--------------'
      a = ( pa * 180. / pi )
c KS CHANGED THIS:
c      if( lvels )a = a - 90.
      a = a - 90.
      write( 4, '( a, 8x, f8.2 )' )'disk PA, phi_d^prime (deg):', a
      write( 4, '( a, 29x, f5.2 )' )'disk eps:', eps
      write( 4, '( a, 11x, 2f8.2 )' )'x,y center (data units):',
     +                                                        xcen, ycen
      if( lnax )then
        if( lvels ) then
          write( 4, '( a, 12x, f8.2 )' )
     +                   'Non-axisymm phib (deg):', ( phib * 180. / pi )
          write( 4, '( a, 16x, i7 )' )'Harmonic order m:', order
        end if
        if( lphot ) then
          write( 4, '( a, 22x, f8.2 )' )
     +                   'Bar PA (deg):', ( bar_pa * 180. / pi - 90. )
          write( 4, '( a, 27x, f8.2 )' )'Bar eps:', bar_eps
        end if
      else
        write( 4, * )
        write( 4, * )
      end if
c Inputs specific to velocities:
      if( lvels )then
        write( 4, '( a, 23x, f8.2 )' )'Vsys (km/s):', vsys
        write( 4, '( a, 18x, f8.2 )' )'Delta_ISM (km/s):', eISM
        if( evelfile )then
          write( 4, '( a, 16x, f8.2 )' )'Delta_D^max (km/s):', errtol
        else
          write( 4, * )
        end if
        if( lwarp )then
          write( 4, '( a, 18x, f8.2 )' )'r_w (data units):', rwarp
          write( 4, '( a, 21x, f8.2 )' )'Warp eps welm:', wepsm
          write( 4, '( a, 21x, f8.2 )' )'Warp PA wphim:', wphim
        else
          write( 4, * )
          write( 4, * )
          write( 4, * )
        end if
      end if
c input values specific to photometry
      if( lphot )then
        if( lbulge )then
          write( 4, '( a, 16x, f8.2 )' )
     +                        'Bulge r_e (pixels):', r_bulge
          write( 4, '( a, 26x, f8.2 )' )'Sersic n:', bulge_n
          write( 4, '( a, 25x, f8.2 )' )'Bulge eps:', bulge_el
        else
          write( 4, * )
          write( 4, * )
          write( 4, * )
        end if
        write( 4, * )
        write( 4, * )
        write( 4, * )
      end if
c write out seeing, smoothing, bootstrap details
      write( 4, * )
      if( lseeing )then
        write( 4, '( a, 2x, f8.2 )' )
     +             'Seeing correction applied. FWHM (data units):', rsee
      else
        write( 4, '( a )' )'No seeing correction applied'
      end if
      if( ( lambda1 .gt. 0. ) .or. ( lambda2 .gt. 0. ) )then
        write( 4, '( a, 2x, f8.2, f8.2 )' )
     +                     'Smoothing params l1, l2:', lambda1, lambda2
      else
        write( 4, '( a )' )'No model component smoothing applied'
      end if
      if( luncert )then
        write( 4, '( a, 2(a7, 2x, i7 ), a7, f8.2 )' )
     +    'Uncertainties estimated via bootstrap:',
     +                 '  seed:',     seed,'  nunc:', nunc,
     +                 '  junc:',     junc
      else
        write( 4, '( a )' ) 'No uncertainties estimated'
      end if
      if( ( lradial .or. lnax ) .and.
     +    ( ( minr .lt. sma0 ) .or. ( maxr .lt. smaf ) ) )then
        write( 4, '( a, 2f8.2 )' )
     +                    'Non-axisymmetries fixed to 0 outside range:',
     +                                        sma( nminr ), sma( nmaxr )
      else
        write( 4, * )
      end if
      write( 4, * )
      write( 4, * )
c output best fitting values of the parameters
      write( 4, '( a )' )'Best fitting values'
      write( 4, '( 21a )' )( dash, i = 1, 21 )
      if( ldisk )then
        if( lpa )then
          astpa = newpa * 180. / pi - 90.
          dmodpa = mod( astpa, 360. )
          write( 4, '( a, 8x, f8.2, a5, f5.2 )' )
     +     'disk PA, phi_d^prime (deg):', dmodpa, pm, enewpa * 180. / pi
        else
          write( 4, * )
          newpa = pa
        end if
        if( leps )then
          newi = acos( 1. - neweps )
          write( 4, '( a, 26x, f8.2, a5, f5.2 )' )'disk eps:',
     +                                             neweps, pm, eneweps
          write( 4, '( a, 19x, f8.2, a5, f5.2 )' )'disk incl (deg):',
     +                                             newi * 180. / pi, pm,
     +                                 eneweps / sin( newi ) * 180. / pi
        else
          write( 4, * )
          neweps = eps
        end if
        if( lcentre )then
          write( 4, '( a, 11x, f8.2, a5, f5.2, a1, f8.2, a5, f5.2 )' )
     +                                      'x,y center (data units):',
     +               newxcen, pm, enewxcen, comma, newycen, pm, enewycen
        else
          write( 4, * )
        end if
      end if
c output parameters for velocities
      if( lvels )then
        if( lnax .and. lphib )then
          naxphi = newphib
          phibsky1 = atan( tan( naxphi )* ( 1. - neweps ) ) + newpa
          phibsky2 =
     +           atan( tan( naxphi + pi/2. ) * ( 1. - neweps ) ) + newpa
          phibsky1 = ( phibsky1 * 180. / pi ) - 90.
          phibsky2 = ( phibsky2 * 180. / pi ) - 90.
          write( 4, '( a, 1x, f8.2, a5, f5.2, a1, f8.2, a1, f8.2 )' )
     +    'Non-axisymm phib (disk plane, deg):', ( naxphi * 180. / pi ),
     +    pm, ( enewphib * 180. / pi ), comma, phibsky1, comma, phibsky2
        else
          write( 4, * )
        end if
        if( lsystemic ) then
          write( 4, '( a, 23x, f8.2, a5, f5.2 )' )
     +                                 'Vsys (km/s):', ivsys, pm, eivsys
        else
          write( 4, * )
        end if
        if( lwarp )then
          if( lrwarp )then
            write( 4, '( a, 18x, f8.2, a5, f5.2 )' )
     +                            'r_w (data units):', newrw, pm, enewrw
          else
            write( 4, * )
          end if
          if( lwepsm )then
            write( 4, '( a, 21x, f8.2, a5, f5.2 )' ) 'Warp eps welm:',
     +                 ( 180. * newpm / pi ), pm, ( 180. * enewpm / pi )
          else
            write( 4, * )
          end if
          if( lwpm )then
            write( 4, '( a, 21x, f8.2, a5, f5.2 )' )
     +                               'Warp PA wphim:', newem, pm, enewem
          else
            write( 4, * )
          end if
        else
          write( 4, * )
          write( 4, * )
          write( 4, * )
        end if
        write( 4, * )
        write( 4, * )
        write( 4, * )
        write( 4, * )
        write( 4, * )
c        write( 4, * )
      end if
c output params, photometry
      if( lphot )then
        if( lnax )then
          if( lphib )then
            astpa = 180. * newbar_pa / pi - 90.
            bmodpa = amod( astpa, 360. )
            ratio = max(1. - neweps, 0.01)
            planepa = atan( tan( newbar_pa - newpa ) / ratio ) 
            planeastpa = ( 180. * planepa / pi ) - 90.
            write( 4, '( a, 1x, f8.2, a5, f5.2, a1, f8.2 )' )
     +             'Non-axisymm phib (sky plane, deg):',
     +             bmodpa, pm, 180. * enewbar_pa / pi, comma, planeastpa
          else
            write( 4, * )
          end if
          if( lepsb )then
            write( 4, '( a, 27x, f8.2, a5, f5.2 )' )
     +                           'Bar eps:', newbar_eps, pm, enewbar_eps
          else
            write( 4, * )
          end if
        else
          write( 4, * )
          write( 4, * )
        end if
        if( lbulge )then
          if( lr_e )then
            write( 4, '( a, 16x, f8.2, a5, f5.2 )' )
     +            'Bulge r_e (pixels):', newr_bulge, pm, enewr_bulge
          else
            write( 4, * )
          end if
          if( lI_0 )then
            write( 4, '( a, 19x, f8.2, a5, f5.2 )' )
     +            'Bulge I_0 (ADU):', ibulge( 1 ), pm, eibulge( 1 )
          else
            write( 4, * )
          end if
          if( lsersn )then
            write( 4, '( a, 26x, f8.2, a5, f5.2 )' )
     +            'Sersic n:',  newbulge_n, pm, enewbulge_n
          else
            write( 4, * )
          end if
          if( lbleps )then
            write( 4, '( a, 25x, f8.2, a5, f5.2 )' )
     +                         'Bulge eps:', newbulge_l, pm, enewbulge_l
          else
            write( 4, * )
          end if
        else
          write( 4, * )
          write( 4, * )
          write( 4, * )
          write( 4, * )
        end if
        write( 4, * )
        if( lnax .or. lbulge )then
          write( 4, '( a, 11x, f8.2, a5, f5.2 )' )
     +   'Disk light fraction (%):', 100. * dskfrac, pm, 100. * edskfrac
          if( lnax )then
            write( 4, '( a, 12x, f8.2, a5, f5.2 )' )
     +    'Bar light fraction (%):', 100. * barfrac, pm, 100. * ebarfrac
          else
            write( 4, * )
          endif
          if( lbulge )then
            write( 4, '( a, 10x, f8.2, a5, f5.2 )' )
     +  'Bulge light fraction (%):', 100. * blgfrac, pm, 100. * eblgfrac
          else
            write( 4, * )
          end if
        else
          write( 4, * )
          write( 4, * )
          write( 4, * )
        end if
      end if
c output minimization details
      write( 4, * )
      write( 4, '( a )' )'Minimization Details '
      write( 4, '( 21a )' )( dash, i = 1, 21 )
      write( 4, '( a, 12x, i7 )' )'# points Dn used in fit:', pixct
      write( 4, '( a, 7x, i7 )' )
     +                       '# iterations in minimization:', final_iter
      write( 4, '( a, 11x, f12.6 )' )'Minimum chi^2 found:', final_chi2
      write( 4, '( a, 10x, i7 )' )
     +                          'Degrees of freedom in fit:', int( dof )
      write( 4, * )
      if( lvels )then
c output fitted velocities
        write( 4, '( a )' )'Fitted velocity components
     +(radii in data units, velocities in km/s):'
        write( 4, '( a )' )'       r        npts          Vt' //
     +                '         eVt          Vr         eVr' //
     +                '        Vm,t       eVm,t        Vm,r       eVm,r'
        write( 4, '( 116a )' )( dash, i = 1, 116 )
        do i = 1, nellip
          write( 4, '( 2f10.2, 8f12.2 )' )sma( i ), ringpts( i ),
     +            id( i ), eid( i ), irad( i ), eirad( i ), ibitan( i ),
     +            eibitan( i ), ibirad( i ), eibirad( i )
        end do
      end if
      if( lphot )then
c output fitted velocities
        write( 4, '( a )' )'Fitted intensities
     +               (radii in pixels, intensities in ADU):'
        write( 4, '( a )' )'       r        npts       Idisk' //
     +                '      eIdisk        Ibar       eIbar' //
     +                '     Ibulge      eIbulge'
        write( 4, '( 116a )' )( dash, i = 1, 93 )
        do i = 1, nellip
          write( 4, '( 2f10.2, 6f12.2 )' )sma( i ), ringpts( i ),
     +           id( i ), eid( i ), ibar( i ), eibar( i ),
     +           blgprof( i ), eblgprof( i )
        end do
      end if
c
      close( 4 )
c report success
      i = lnblnk( outpfile )
      if( ltext )then
        j = lnblnk( outroot )
        print *,
     +       'Wrote out: ', outpfile( 1:i ), '  ', outroot( 1:j )//'mod'
      else
        print *, 'Wrote out: ', outpfile( 1:i )
      end if
      return
      end
