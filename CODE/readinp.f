      subroutine readinp
c This subroutine reads the input file. The structure of the input file
c   is given in the documentation.
c
c   written by KS
c   Adapted by RZS summer 2009
c   Polished by JAS October 2009
c   Incorporated photometric fitting by JAS March 2011
c   Simple warp option added by JAS July 2011
c   Polished KS March 2012
c   Added call to initlz - JAS June 2012
c
c   lines of characters in the input file can be either in single quotes
c   or plain text.  If the first character of a line is ' [char(39)], the
c   string is presumed to start at the second character
c
c common blocks
      include 'commons.h'
c
c externals
      integer lnblnk, neari
c
c local variables:
      integer i, ic, istat, j, k
      character line*120, str*3
c
c we always fit for a disk or the circular speed
      ldisk = .true.
c initialize other logical variables
      leps      = .false.
      lcentre   = .false.
      lpa       = .false.
      lnax      = .false.
      lphib     = .false.
      lepsb     = .false.
      lbulge    = .false.
      lbleps    = .false.
      lI_0      = .false.
      lr_e      = .false.
      lsersn    = .false.
      lwarp     = .false.
      lrwarp    = .false.
      lwepsm    = .false.
      lwpm      = .false.
      linter0   = .false.
      lradial   = .false.
      lseeing   = .false.
      lsystemic = .false.
      luncert   = .false.
      lFITS     = .false.
      ltext     = .false.
      VELmps    = .false.
c try to open specified file
      open( unit = 13, file = infile, status = 'old', iostat = i )
      if( i .ne. 0 )then
        print '( a, 2x, a )', 'Failed to find input file:', infile
        call crash( 'readinp' )
      end if
c skip a line
      read( 13, * )
c initialize line counter
      ic = 1
c l2: specify whether the fit is to a photometric image or to a velocity map
      call getline( 13, line, istat, ic )
      if( line( 1:1 ) .eq. char( 39 ) )then
        str = line( 2:4 )
      else
        str = line( 1:3 )
      end if
c convert to upper case if not already
      do i = 1, 3
        j = ichar( str( i:i ) )
        if( ( j .gt. 96 ) .and. ( j .lt. 123 ) )then
          str( i:i ) = char( j - 32 )
        end if
      end do
      lvels = str .eq. 'VEL'
      lphot = str .eq. 'PHO'
      if( .not. ( lvels .or. lphot ) )go to 1
c l3: FITS type data flag - line ignored for lphot
      call getline( 13, line, istat, ic )
      if( lvels )then
        read( line, *, err = 1, end = 1 )lFITS, VELmps
        lText = .not. lFITS
      else
c photometry data are assumed in FITS format
        lFITS = .true.
      end if
c l4: name of file with data
      call getline( 13, line, istat, ic )
      if( lvels )then
        if( line( 1:1 ) .eq. char( 39 ) )then
          read( line, *, err = 1, end = 1 )invfile
        else
          invfile = line
        end if
      else
        if( line( 1:1 ) .eq. char( 39 ) )then
          read( line, *, err = 1, end = 1 )inpfile
        else
          inpfile = line
        end if
      end if
c l5:
      call getline( 13, line, istat, ic )
      if( lvels )then
c name of file that contains velocity uncertainties
        if( line( 1:1 ) .eq. char( 39 ) )then
          read( line, *, err = 1, end = 1 )inevfile
        else
          inevfile = line
        end if
      else
c parameters of photometry
        read( line, *, err = 1, end = 1 )sky, skysig, gain
      end if
c l6: read if lFITS - otherwise ignored
      call getline( 13, line, istat, ic )
      if( lFITS )then
c values for selecting the region of the image to fit
        read( line, *, err = 1, end = 1 )lox, loy, xrange, yrange
        inp_pts = xrange * yrange
      end if
c l7: FITS sampling parameters - ignored if .not. lFITS
      call getline( 13, line, istat, ic )
      if( lFITS )then
c read image region parameters needed only for FITS input
        read( line, *, err = 1, end = 1 )regrad, regpa, regeps,
     +        istepout, pixscale
      else
        pixscale = 1
      end if
c l8: name of file for parameter output
      call getline( 13, line, istat, ic )
      if( line( 1:1 ) .eq. char( 39 ) )then
        read( line, *, err = 1, end = 1 )outpfile
      else
        outpfile = line
      end if
c parse name of output file for subdirectories and root name
      i = index( outpfile, '/' )
      j = lnblnk( outpfile )
      k = i + 1
      do while ( i .gt. 0 )
        i = index( outpfile( k:j ), '/' )
        k = k + i
      end do
      i = index( outpfile( k:j ), '.' )
      k = k + i - 1
      outroot = outpfile( 1:k )
c l9: toggles for axisymmetric part of fit
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )lpa, leps, lcentre
c l10: initial guesses for disk PA and eps = (1 - b/a)
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )pa, eps
c l11: initial guess for galaxy center
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )xcen, ycen
c l12: toggles for non-axisymmetric fitting
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )lnax
      if( lnax )then
        if( lvels )then
          read( line, *, err = 1, end = 1 )lnax, lphib, phib, order
c check for an allowed value of order
          if( order .lt. 1 .or. order .gt. 2 )then
            print *, 'Order must be 1 or 2'
            go to 1
          end if
        else
c bar position angle and apparent ellipticity (1-b/a)
          read( line, *, err = 1, end = 1 )lnax, lphib, lepsb,
     +                                     bar_pa, bar_eps
        end if
      else
c supply defaults
        lphib = .false.
        phib = 0
        if( lvels )then
          order = 2
        else
          lepsb = .false.
          bar_pa = 0
          bar_eps = 0.5
        end if
      end if
c l13: ignored if phot
      call getline( 13, line, istat, ic )
      if( lvels )then
        read( line, *, err = 1, end = 1 )linter0, lradial
      end if
c l14:
      call getline( 13, line, istat, ic )
      if( lvels )then
c initial systemic velocity, ISM turbulence parameter, & vel error tolerance
        read( line, *, err = 1, end = 1 )lsystemic, vsys, eISM, errtol
      else
        read( line, *, err = 1, end = 1 )lbulge
        if( lbulge )then
c some individual bulge parameter toggles and values
          read( line, *, err = 1, end = 1 )lbulge, lr_e, r_bulge
c always fit for I_0 for now
          lI_0 = .true.
        else
c supply defaults
          lr_e = .false.
          rbulge = 1
          lI_0 = .false.
          bulge_I0 = 0
        end if
      end if
c l15:
      call getline( 13, line, istat, ic )
      if( lvels )then
c warp toggle and parameters
        read( line, *, err = 1, end = 1 )lwarp
c set invidual warp toggles and initial parameters
        if( lwarp )then
          read( line, *, err = 1, end = 1 )lwarp, lrwarp, lwepsm, lwpm,
     +                                     rwarp, wepsm, wphim
        else
c ignore rest of the line if lwarp is false - values are unused
          lrwarp = .false.
          lwepsm = .false.
          lwpm = .false.
          rwarp = 1
          wphim = 0
          wepsm = 0
        end if
      else
c more bulge toggles and initial values
        if( lbulge )then
          read( line, *, err = 1, end = 1 )lsersn, lbleps,
     +                                     bulge_n, bulge_el
        else
c set defaults
          lsersn = .false.
          lbleps = .false.
          bulge_n = 1
          bulge_el = 0
        end if
      end if
c l16: seeing correction - rsee = FWHM of seeing disk
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )rsee
      lseeing = rsee .gt. 0.
c l17: smoothing parameters
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )lambda1, lambda2
c l18: uncertainty parameters
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )luncert
c read bootstrap parameters only if uncertainties were requested
      if( luncert )read( line, *, err = 1, end = 1 )luncert,
     +                                                  seed, nunc, junc
c l19: verbose toggle
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )verbose
c l20:
      call getline( 13, line, istat, ic )
      if( lnax .or. lradial )then
c minimum, maximum radius for bar or non-circular flow fit
        read( line, *, err = 1, end = 1 )minr, maxr
      else
c supply defaults
        minr = 0
        maxr = 10000
      end if
c l21: input ring radii
      istat = 0
      i = 0
      do while ( istat .eq. 0 )
        i = i + 1
        ic = ic + 1
        read( 13, *, iostat = istat )sma( i )
      end do
      close( 13 )
c finished reading file
      nellip = i - 1
      sma0 = sma( 1 )
      smaf = sma( nellip )
c may need to adjust nminr & nmaxr
      if( ( maxr .lt. smaf ) .and. ( maxr .ge. sma0 ) )then
        i = neari( sma, nellip, maxr )
        maxr = sma( i )
        nmaxr = i
      else if( maxr .ge. smaf )then
        nmaxr = nellip
        maxr = smaf
      else
        nmaxr = 0
      end if
      if( ( minr .lt. maxr ) .and. ( minr .ge. sma0 ) )then
        i = neari( sma, nellip, minr )
        minr = sma( i )
        nminr = i
      else if( minr .le. sma0 )then
        nminr = 1
        minr = sma0
      end if
c zero or negative radial range forces a simple model with only circular flow
      if( nminr .gt. nmaxr )then
        lradial = .false.
        lnax = .false.
      end if
c set up parameters for seeing correction if selected
      if( lseeing )call blurset
c
c compute inclination from input ellipticity: cos( incl ) = 1 - eps
      incl = acos( 1. - eps ) * 180. / pi
c change to definition PA used in program (y-axis to left)
      pa = pa + 90.
      if( lphot .and. lnax )bar_pa = bar_pa + 90.
c convert to radians
      pa = pa * pi / 180.
      bar_pa = bar_pa * pi / 180.
      phib = phib * pi / 180.
      wphim = wphim * pi / 180.
c set up nd and ntot
      call initlz
      return
c read failed
    1 print '( 2x, a )', line
      print *, 'Problem with line', ic
      print '( ''Error reading file '', a )', infile
      call crash( 'readinp' )
      end
