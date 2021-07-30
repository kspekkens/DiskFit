c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c Main driving program for fitting data
c
c This software package was developed over a period of years from 2002 by
c   Jerry Sellwood and Kristine Spekkens in collaboration with:
c   Eric I Barnes, Adam Reese, and Ricardo Zanmar Sanchez
c
c Fits velocity fields of galaxies.   - KS June 2008
c Add extra flag and allow initial biphib
c Split disk flag into lPA & leps     - RZS April 2009
c Polished and smoothing option added - JAS October 2009
c AR photometry fitting option resurrected and
c     seeing correction implemented   - JAS March 2011
c Warp model added                    - JAS July 2011
c Released as Diskfit 1.0             - Oct 2011
c upgraded to f90                     - Jan 2015
c Fixed bugs:
c   in seeing corrections             - Jul 2015
c   in rastering of text input        - Aug 2015
c   output selected part of an image  - Sep 2015
c Phot allow user input uncertainties - Sep 2015
c Released as DiskFit 1.2             - Sep 2015
c
      include 'commons.h'
c
c local variables
      integer i, j
      real*8 chi2, eparams( md ), params( md )
c
      print *, 'DiskFit (v1.2.2)  Copyright (C) 2015, 2017,' //
     +         ' Jerry Sellwood and Kristine Spekkens'
      print *
      print *, 'This program comes with ABSOLUTELY NO WARRANTY.' //
     +         '  It is free software'
      print *, 'and you are welcome to' //
     +         ' redistribute it under certain conditions.'
      print *, 'See <http://www.gnu.org/licenses/> for details.'
      print *
c read input values
      print *,'Enter input parameter file name'
      read( *, * )infile
      call readinp
      if( lwarp .and. lphot )then
        print *, 'Warp model for photometry not programmed'
        call crash( 'Main' )
      end if
c initialize error vectors
      do i = 1, nd
        eparams( i ) = 0
      end do
      do i = 1, ntot
        efitval( i ) = 0
      end do
c read in image or map
      call prepdata
      if( l2d )then
        i = xrange
        j = yrange
      else
        i = inp_pts
        j = 1
      end if
c
c set initial guesses for non-linear parameters
      call setpars( params )
c
c screen output before minimization
      call screenoutput
c allocate space
      if( l2d )then
        allocate ( iwgt( mk, xrange, yrange ) )
        allocate ( wgt( mk, xrange, yrange ) )
      else
        allocate ( iwgt( mk, inp_pts, 1 ) )
        allocate ( wgt( mk, inp_pts, 1 ) )
      end if
c
c find minimum of chi^2
      chi2 = 0.d0
      iter = 0
      call mini( params, chi2, iter )
      final_chi2 = chi2
      final_iter = iter
c
c evaluate best fitting solution
      call bestfit( params, eparams )
c create model and residuals images or maps - residuals need for bootstrap
      call bestmod( params )
c output FITS files of best-fit model with all pixel values in working box
      if( lFITS )call writemod
c
c run a bootstrap loop if desired
      if( luncert )call bootstrap( params, eparams )
c
c write out optimal parameters and values, with their estimated uncertainties
      call writeout
      end
