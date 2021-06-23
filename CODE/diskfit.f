c Main driving program for fitting data
c
c This software package was developed over a period of years from 2002 by
c   Jerry Sellwood in collaboration with:
c   Eric I Barnes,
c   Adam Reese,
c   Kristine Spekkens, and
c   Ricardo Zanmar Sanchez
c
c Fits velocity fields of galaxies.   - KS June 2008
c Add extra flag and allow initial biphib
c Split disk flag into lPA & leps     - RZS April 2009
c Polished and smoothing option added - JAS October 2009
c AR photometry fitting option resurrected and
c     seeing correction implemented   - JAS March 2011
c Warp model added                    - JAS July 2011
c Released as Diskfit 1.0             - Oct 2011
c
      include 'commons.h'
c local variables
      integer i
      real*8 chi2, eparams( md ), params( md )
c
c read input values
      print *
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
c
c set initial guesses for non-linear parameters
      call setpars( params )
c
c screen output before minimization
      call screenoutput
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
