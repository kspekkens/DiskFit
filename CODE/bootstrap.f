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
c
      include 'commons.h'
c
c calling arguments
      real*8 params( md ), eparams( md )
c
c external
      integer lnblnk
c
c local arrays
      real, allocatable :: sdot( :, : )
      real, allocatable :: x2ros( :, : )
      real*8, allocatable :: afitval( :, : )
      real*8, allocatable :: aparams( :, : )
      real*8, allocatable :: inparams( : )
      real*8, allocatable :: fitval0( : )
      real*8, allocatable :: lfrac( :, : )
      real*8, allocatable :: params0( : )
c
c local variables
      character form*7
      integer i, ierr, ii, ix, iy, j, jdum, k, ntimes
      real barfrac_loc, blgfrac_loc, dskfrac_loc
      real*8 chi2
      character*100 erasefile, bootfile
c
      allocate ( afitval( nunc, ntot ) )
      allocate ( aparams( nunc, ntot ) )
      allocate ( inparams( nd ) )
      allocate ( fitval0( ntot ) )
      allocate ( lfrac( 3, nunc ) )
      allocate ( params0( nd ) )
c set up blurring arrays in case they were not preset - no harm to do it again
      if( l2D )call blurset
c
      jdum = seed
c
c save best fitting parameters, initial data, and best fit x2res map
      do i = 1, nd
        params0( i ) = params( i )
      end do
      do j = 1, ntot
        fitval0( j ) = fitval( j )
      end do
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
c save best fitting light fractions
      if( lphot )then
        dskfrac_loc = dskfrac
        barfrac_loc = barfrac
        blgfrac_loc = blgfrac
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
c save results to output file for post-proecessing
      open( 4, file = bootfile, status = 'unknown' )
      write( 4, '( a )' )infile
      write( 4, '( 7( a9, l2 ) )' )' lpa: ', lpa, ' leps:', leps,
     *             ' lcentre: ', lcentre, ' lsystemic:   ', lsystemic,
     *             ' lradial: ', lradial, ' lnax:        ', lnax,
     *             ' luncert: ', luncert
      write( 4, '( a )' )
     +               '#nunc   npars    ntot   nellip  nminr  nmaxr seed'
      write( 4, '( 7i7 )' )nunc, nd, ntot, nellip, nminr, nmaxr, seed
      write( 4, '( a )') ' '
c output best-fit parameters to bootstrap file
      ii = 0
      write( 4, * )'Bootstrap number: ', ii
      write( 4, * )( params( j ), j = 1, nd )
      write( 4, * )( fitval( j ), j = 1, ntot )
      if( lphot .and. ldisk .and. ( lnax .or. lbulge ) )
     +                            write( 4, * )dskfrac, barfrac, blgfrac
      write( 4, * )
c create flag file - a graceful exit occurs should this file be deleted
      open( 1, file = erasefile, status = 'UNKNOWN' )
      write( 1, '( a )') 'Erase this file to exit DiskFit gracefully'
      close( 1 )
c bootstrap loop
      ierr = 0
      ii = 0
      do while ( ierr .eq. 0 )
        ii = ii + 1
        print *
        print *
        print *, 'Bootstrap number: ', ii
c re-initialize params
        do j = 1, nd
          params( j ) = inparams( j )
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
c generate new velocity field - with coherent patches of residuals as SS07
          if( l2D )then
            call pseudp( jdum )
          else
            call pseudo( jdum )
          end if
        else
c generate new velocity field - with rotation and radial scaling as SZS09
          call pseudr( jdum )
        end if
c run minimization
        chi2 = 0.d0
        iter = 0
        call mini( params, chi2, iter )
        if( lnax )then
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
          aparams( ii, j ) = params( j )
        end do
        do j = 1, ntot
          afitval( ii, j ) = fitval( j )
        end do
        if( lphot )then
          lfrac( 1, ii ) = dskfrac
          lfrac( 2, ii ) = barfrac
          lfrac( 3, ii ) = blgfrac
        end if
c output new fitted parameters to bootstrap file
        write( 4, * )'Bootstrap number: ', ii
        write( 4, * )( params( j ), j = 1, nd )
        write( 4, * )( fitval( j ), j = 1, ntot )
        if( lphot .and. ldisk .and. ( lnax .or. lbulge ) )
     +      write( 4, * )lfrac( 1, ii ), lfrac( 2, ii ), lfrac( 3, ii )
        write( 4, * )
c try to open flag file - error on opening is the signal to end gracefully
        open( 3, file = erasefile, status = 'old', iostat = ierr )
        if( ierr .eq. 0 )read( 3, *, iostat = ierr )
        close( 3 )
c end bootstrap loop here
        if( ii .eq. nunc )ierr = 1
      end do
c close output file
      close( 4 )
c update the actual number of bootstrap iterations
      nunc = ii
      print *
      print *
      print *, 'Finished', ii, ' bootstrap iterations'
c restore best fitting parameters and values
      do i = 1, nd
        params( i ) = params0( i )
      end do
      do j = 1, ntot
        fitval( j ) = fitval0( j )
      end do
      if( lphot )then
        dskfrac = dskfrac_loc
        barfrac = barfrac_loc
        blgfrac = blgfrac_loc
      end if
c get error estimates
      call geterrs( aparams, afitval, lfrac, eparams, nunc )
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
