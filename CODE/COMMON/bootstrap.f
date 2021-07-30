      subroutine bootstrap( params, eparams )
c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c estimates uncertainties for best-fit parameters using a bootstrap technique
c  Originally written by KS
c  Modified by RZS in July 2009
c  Polished by JAS October 2009
c  Modified to include photometry by JAS October 2011
c  Fixed problem of angle errors and split off geterrs - JAS June 2012
c  Updated to f90 - JAS Jan 2015
c  Fixed array dimensions bug passed to geterrs,
c     changed to adjust the fitted bar angle for vels only,
c     added calls to fix_angles and tabcov, and
c     close and reopen the .bstrp file every 10 bootstraps - JAS Aug 2017
c  Reordered subscripts of arrays afitval and aparams - JAS Dec 2017
c  Include best fit params in error estimate and csv file - JAS Aug 2020
c  Initialize warp deprojection, if needed - JAS Sep 2020
c
      include 'commons.h'
c
c calling arguments
      real*8 params( md ), eparams( md )
c
c external
      integer lnblnk
c$$$      real*8 ran1_dbl
c
c local arrays
      real, allocatable :: sdot( :, : )
      real, allocatable :: x2ros( :, : )
      real, allocatable :: wparr( :, : )
      real*8, allocatable :: afitval( :, : )
      real*8, allocatable :: aparams( :, : )
      real*8, allocatable :: inparams( : )
      real*8, allocatable :: lfrac( :, : )
c
c local variables
      character form*7
      integer i, ierr, ii, ix, iy, j, jdum, k, ntimes
      real*8 chi2
      character*100 erasefile, bootfile
c
c the first dimension of aparams has to be the same as that of afitval
      allocate ( afitval( ntot, nunc + 1 ) )
      allocate ( aparams( ntot, nunc + 1 ) )
      allocate ( inparams( nd ) )
      allocate ( lfrac( 3, nunc + 1 ) )
c set up blurring arrays in case they were not preset - no harm to do it again
      if( l2D )call blurset
c
      jdum = seed
c
c save best fitting parameters, initial data, and best fit x2res map
      do j = 1, nd
        aparams( j, 1 ) = params( j )
      end do
      do j = 1, ntot
        afitval( j, 1 ) = fitval( j )
      end do
      if( lphot )then
        lfrac( 1, 1 ) = dskfrac
        lfrac( 2, 1 ) = barfrac
        lfrac( 3, 1 ) = blgfrac
      end if
      if( l2D )then
        allocate ( sdot( xrange, yrange ) )
        allocate ( x2ros( xrange, yrange ) )
        do iy = 1, yrange
          do ix = 1, xrange
            sdot( ix, iy ) = sdat( ix, iy )
            x2ros( ix, iy ) = x2res( ix, iy )
          end do
        end do
      else
        allocate ( sdot( inp_pts, 1 ) )
        allocate ( x2ros( inp_pts, 1 ) )
        do i = 1, inp_pts
          sdot( i, 1 ) = sdat( i, 1 )
          x2ros( i, 1 ) = x2res( i, 1 )
        end do
      end if
c set initial guesses for non-linear parameters and keep a copy
      call setpars( params )
      do j = 1, nd
        inparams( j ) = params( j )
      end do
c
      print *
      print *, 'Determining uncertainties...'
      print *, 'Beginning bootstrap'
c construct names for work files
      k = lnblnk( outroot ) + 1
      erasefile = outroot
      write( erasefile( k:k+4 ), '( a5 )' )'erase'
      bootfile = outroot
      write( bootfile( k:k+4 ), '( a5 )' )'bstrp'
c ensure a different file name for each seed in case of parallel bootstrap runs
      i = log10( real( abs( seed ) ) + .5 ) + 1
      i = max( i, 1 )
      form = '( i   )'
      if( i .ge. 10 )then
        write( form( 4:5 ), '( i2 )' )i
      else
        write( form( 4:4 ), '( i1 )' )i
      end if
      k = k + 5
      write( bootfile( k:k+i-1 ), form )abs( seed )
c save results to a new output file for post-processing
      open( 4, file = bootfile, status = 'new', iostat = ii )
      if( ii .ne. 0 )then
        print *, 'failed to open new bootstrap file'
        call crash( 'BOOTSTRAP' )
      end if
      write( 4, '( a )' )infile
      write( 4, '( 7( a9, l2 ) )' )' lpa: ', lpa, ' leps:', leps,
     *             ' lcentre: ', lcentre, ' lsystemic:   ', lsystemic,
     *             ' lradial: ', lradial, ' lnax:        ', lnax,
     *             ' luncert: ', luncert
      write( 4, '( a )' )
     +               '#nunc   npars    ntot   nellip  nminr  nmaxr seed'
      write( 4, '( 7i7 )' )nunc, nd, ntot, nellip, nminr, nmaxr, seed
      write( 4, * )
c output best-fit parameters to bootstrap file
      ii = 0
      write( 4, * )'Bootstrap number: ', ii, ' (Best fit)'
      write( 4, * )( aparams( j, 1 ), j = 1, nd )
      write( 4, * )( sngl( afitval( j, 1 ) ), j = 1, ntot )
      if( lphot .and. ldisk .and. ( lnax .or. lbulge ) )
     +                            write( 4, * )dskfrac, barfrac, blgfrac
      write( 4, * )
c create flag file - a graceful exit occurs should this file be deleted
      open( 1, file = erasefile, status = 'unknown' )
      write( 1, '( a )') 'Erase this file to exit DiskFit gracefully'
      close( 1 )
c setup and save warp arrays, if needed
      if( lwarp )then
        call warpinit( params )
        allocate ( wparr( 6, nellip ) )
        do i = 1, nellip
          wparr( 1, i ) = wel( i )
          wparr( 2, i ) = wphi( i )
          wparr( 3, i ) = wba( i )
          wparr( 4, i ) = wcp( i )
          wparr( 5, i ) = wsp( i )
          wparr( 6, i ) = wsi( i )
        end do
      end if
      if( nint( junc ) .ge. mblur )print '( a, i5 )', 'Warning: proper'
     + // ' functioning of this bootstrap option requires junc <', mblur
c bootstrap loop
      ierr = 0
      ii = 1
      do while ( ierr .eq. 0 )
        ii = ii + 1
        print *
        print *
        print *, 'Bootstrap number: ', ii - 1
c re-initialize params
        do j = 1, nd
          params( j ) = inparams( j )
c$$$          params( j ) = inparams( j ) *
c$$$     +                           ( 1 + .04 * ( ran1_dbl( jdum ) - .5 ) )
        end do
c restore input x2res map
        if( l2D )then
          do iy = 1, yrange
            do ix = 1, xrange
              x2res( ix, iy ) = x2ros( ix, iy )
            end do
          end do
        else
          do i = 1, inp_pts
            x2res( i, 1 ) = x2ros( i, 1 )
          end do
        end if
        if( junc .ge. 0.0 ) then
c generate new map - with coherent patches of residuals as SS07
          if( l2D )then
            call pseudp( jdum )
          else
            call pseudo( jdum )
          end if
        else
c restore warp arrays, if needed
          if( lwarp )then
            do i = 1, nellip
              wel( i ) = wparr( 1, i )
              wphi( i ) = wparr( 2, i )
              wba( i ) = wparr( 3, i )
              wcp( i ) = wparr( 4, i )
              wsp( i ) = wparr( 5, i )
              wsi( i ) = wparr( 6, i )
            end do
          end if
c generate new map - with rotation and radial scaling as SZS09
          call pseudr( jdum )
        end if
c run minimization
        chi2 = 0.d0
        iter = 0
        call mini( params, chi2, iter )
c allow for fit having returned the angle of the minor axis - vels only
        if( lvels .and. lnax )then
          j = 0
          if( lpa )j = j + 1
          if( leps )j = j + 1
          if( lcentre )j = j + 2
          j = j + 1
          k = 0
          k = k + nellip
          if( lradial )k = k + nradial
          ntimes = int( params( j ) * 2.d0 / pi )
          if( ntimes .gt. 0 )then
            params( j ) = params( j ) - dble( ntimes ) * pi / 2.d0
            do j = 1, nnasymm
              fitval( k + j ) =
     +              ( -1.d0 )**abs( ntimes ) * fitval( k + j )
              fitval( k + nnasymm + j ) =
     +              ( -1.d0 )**abs( ntimes ) * fitval( k + nnasymm + j )
            end do
          end if
        end if
c save outputs in arrays
        do j = 1, nd
          aparams( j, ii ) = params( j )
        end do
        do j = 1, ntot
          afitval( j, ii ) = fitval( j )
        end do
        if( lphot )then
          lfrac( 1, ii ) = dskfrac
          lfrac( 2, ii ) = barfrac
          lfrac( 3, ii ) = blgfrac
        end if
c output new fitted parameters to bootstrap file
        write( 4,
     + '( '' Bootstrap number: '', i5, ''/'', i5, ''  Chisq '', f12.5 )'
     +                                               )ii - 1, nunc, chi2
        write( 4, * )( params( j ), j = 1, nd )
        write( 4, * )( fitval( j ), j = 1, ntot )
        if( lphot .and. ldisk .and. ( lnax .or. lbulge ) )
     +      write( 4, * )lfrac( 1, ii ), lfrac( 2, ii ), lfrac( 3, ii )
        write( 4, * )
c flush output buffer after every 10 bootstraps
        if( mod( ii, 10 ) .eq. 0 )then
          close( 4 )
          open( 4, file = bootfile, status = 'old', access = 'append' )
        end if
c try to open flag file - error on opening is the signal to end gracefully
        open( 3, file = erasefile, status = 'old', iostat = ierr )
        if( ierr .eq. 0 )read( 3, *, iostat = ierr )
        close( 3 )
c end bootstrap loop here
        if( ii .eq. nunc + 1 )ierr = 1
      end do
c close output file and report
      close( 4 )
      print *
      print *
      print *, 'Finished', ii - 1, ' bootstrap iterations'
c restore best fitting parameters and values
      do i = 1, nd
        params( i ) = aparams( i, 1 )
      end do
      do j = 1, ntot
        fitval( j ) = afitval( j, 1 )
      end do
      if( lphot )then
        dskfrac = lfrac( 1, 1 )
        barfrac = lfrac( 2, 1 )
        blgfrac = lfrac( 3, 1 )
      end if
c get error estimates
      i = nunc + 1
      nunc = ii
      if( nunc .gt. 1 )then
        call fix_angles( aparams, ntot, i )
        call tabcov( aparams, afitval, ntot, lfrac, i )
        call geterrs( aparams, afitval, ntot, lfrac, eparams, i )
      end if
c restore intial data
      if( l2D )then
        do iy = 1, yrange
          do ix = 1, xrange
            sdat( ix, iy ) = sdot( ix, iy )
          end do
        end do
      else
        do i = 1, inp_pts
          sdat( i, 1 ) = sdot( i, 1 )
        end do
      end if
c reset the weights for the original parameters
      fronly = .false.
      call setwgt( params, lseeing, k )
c save uncertainty estimates
      call bestfit( params, eparams )
      return
      end
