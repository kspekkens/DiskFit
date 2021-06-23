      subroutine mini( params, chi2, iterl )
c this subroutine implements the Powell chi-squared minimization routine
c   and returns the best fit parameters as well as the minimum chi-squared
c   value chi2
c
c   Provenance unknown
c   Polished by JAS Oct 2009
c
      include 'commons.h'
c
c calling arguments
      integer iterl
      real*8 chi2, params( md )
c
c external
      real*8 func
c
c local variables
      integer i, j
      real*8 ftol, xi( md, md )
c
      iterl = 0
c call func only if no non-linear params
      if( nd .eq. 0 )then
        print *
        chi2 = func( params )
      else
c Powell minimization
        ftol = 1.0d-5
c set direction matrix along parameter axes - Powell will scale appropriately.
        do j = 1, nd
          do i = 1, nd
            xi( i, j ) = 0
            if( i .eq. j )xi( i, j ) = 1
          end do
        end do
c main call
        call powell_dbl( params, xi, nd, md, ftol, iterl, chi2 )
        print *
        print *, 'Done minimization'
      end if
      print *, 'Total number of iterations: ', iterl
      print *, 'Minimum chi^2 found: ', sngl( chi2 )
      print *
c estimate light fractions in each component
      if( lphot )call lfracs( params )
      return
      end
