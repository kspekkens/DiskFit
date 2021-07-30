      subroutine readinp
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
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
c   Modified to allow image mask - JAS June 2013
c   Allow for a file of uncertainties for PHOT - JAS Sep 15
c   Initialized order=2 for a bar for PHOT - JAS Mar 17
c   Set skpbimod as 3rd parameter vels option on line 3 - JAS Aug 2020
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
      character line*120, str*4
c
c we always fit for a disk or the circular speed
      ldisk     = .true.
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
c line 2: specify whether the fit is to a photometric image or to a velocity map
      call getline( 13, line, istat, ic )
      if( line( 1:1 ) .eq. char( 39 ) )then
        str( 1:3 ) = line( 2:4 )
      else
        str( 1:3 ) = line( 1:3 )
      end if
c convert to upper case if not already
      do i = 1, 3
        j = ichar( str( i:i ) )
        if( ( j .gt. 96 ) .and. ( j .lt. 123 ) )then
          str( i:i ) = char( j - 32 )
        end if
      end do
      lvels = str( 1:3 ) .eq. 'VEL'
      lphot = str( 1:3 ) .eq. 'PHO'
      if( .not. ( lvels .or. lphot ) )go to 1
c line 3: FITS type data flag - line ignored for lphot
      call getline( 13, line, istat, ic )
      if( lvels )then
        read( line, *, iostat = i )lFITS, VELmps, skpbimod
        if( i .ne. 0 )then
          read( line, *, err = 1, end = 1 )lFITS, VELmps
          skpbimod = .false.
        end if
        lText = .not. lFITS
      else
c photometry data are assumed in FITS format
        lFITS = .true.
c determine whether photometry data come with a separate mask file
        if( line( 1:1 ) .eq. char( 39 ) )then
          str = line( 2:5 )
        else
          str = line( 1:4 )
        end if
c convert to upper case if not already
        do i = 1, 4
          j = ichar( str( i:i ) )
          if( ( j .gt. 96 ) .and. ( j .lt. 123 ) )then
            str( i:i ) = char( j - 32 )
          end if
        end do
        lmask = ( str .ne. 'NONE' ) .and. ( str .ne. '    ' )
c read name of mask file
        if( lmask )then
          if( line( 1:1 ) .eq. char( 39 ) )then
            read( line, *, err = 1, end = 1 )mskfile
          else
            mskfile = line
          end if
        end if
      end if
c line 4: name of file with data
      call getline( 13, line, istat, ic )
      if( line( 1:1 ) .eq. char( 39 ) )then
        read( line, *, err = 1, end = 1 )datfile
      else
        datfile = line
      end if
c line 5:
      call getline( 13, line, istat, ic )
c name of file that contains uncertainties
      if( line( 1:1 ) .eq. char( 39 ) )then
        read( line, *, err = 1, end = 1 )errfile
      else
        errfile = line
      end if
c no file is flagged by the name "NONE" or "    " as the first 4 characters
      str = errfile
c convert to upper case if not already
      do i = 1, 4
        j = ichar( str( i:i ) )
        if( ( j .gt. 96 ) .and. ( j .lt. 123 ) )then
          str( i:i ) = char( j - 32 )
        end if
      end do
      lerrfile =  ( str .ne. 'NONE' ) .and. ( str .ne. '    ' )
      if( lerrfile )then
        if( skpbimod )then
          print *, 'Pixels flagged as bimodal will be skipped'
        else
          print *,
     +          'Pixels flagged as bimodal will be included in the fit'
        end if
      else
        if( skpbimod )then
          print *, 'No error file, so bimodal profiles are not flagged'
          skpbimod = .false.
        end if
      end if
c line 6: read if lFITS - otherwise ignored
      call getline( 13, line, istat, ic )
      if( lFITS )then
c values for selecting the region of the image to fit
        read( line, *, err = 1, end = 1 )lox, loy, xrange, yrange
        inp_pts = xrange * yrange
      end if
c line 7: FITS sampling parameters - ignored if .not. lFITS
      call getline( 13, line, istat, ic )
      if( lFITS )then
c read image region parameters needed only for FITS input
        read( line, *, err = 1, end = 1 )regrad, regpa, regeps,
     +        istepout, pixscale
      else
        pixscale = 1
      end if
c line 8: name of file for parameter output
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
c line 9: toggles for axisymmetric part of fit
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )lpa, leps, lcentre
c line 10: initial guesses for disk PA and eps = (1 - b/a)
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )pa, eps
c line 11: initial guess for galaxy center
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )xcen, ycen
c line 12: toggles for non-axisymmetric fitting
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
c a bar is bisymmetric
          order = 2
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
c line 13: ignored in photometry, unless there is no uncertainties file
      call getline( 13, line, istat, ic )
      if( lvels )then
        read( line, *, err = 1, end = 1 )linter0, lradial
      else
c parameters of photometry
        if( .not. lerrfile )then
          read( line, *, iostat = i )sky, skysig, gain
          if( i .ne. 0 )then
            print *, 'You specified no uncertainties file, so l13' //
     +  ' must contain the parameters needed to compute the noise level'
            go to 1
          end if
        end if
      end if
c line 14:
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
c line 15:
      call getline( 13, line, istat, ic )
      if( lvels )then
c warp toggle and parameters
        read( line, *, err = 1, end = 1 )lwarp
c set individual warp toggles and initial parameters
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
c line 16: seeing correction - dsee = FWHM of seeing disk
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )dsee
      lseeing = dsee .gt. 0.
c line 17: smoothing parameters
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )lambda1, lambda2
c line 18: uncertainty parameters
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )luncert
c read bootstrap parameters only if uncertainties were requested
      if( luncert )read( line, *, err = 1, end = 1 )luncert,
     +                                                  seed, nunc, junc
c line 19: verbose toggle
      call getline( 13, line, istat, ic )
      read( line, *, err = 1, end = 1 )verbose
c line 20:
      call getline( 13, line, istat, ic )
      if( lnax .or. lradial )then
c minimum, maximum radius for bar or non-circular flow fit
        read( line, *, err = 1, end = 1 )minr, maxr
      else
c supply defaults
        minr = 0
        maxr = 10000
      end if
c line 21: input ring radii
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
      incl = sngl( acos( dble( 1. - eps ) ) * 180. / pi)
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
