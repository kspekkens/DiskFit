      subroutine initlz
c sets   nd: the number of non-linear parameters in the fit
c  and ntot: the rank (number of rows) of the matrix
c
c created by JAS from code that was in DiskFit - June 2012
c
c common block
      include 'commons.h'
c
c compute numbers of linear + non-linear parameters
      nd = 0
      ntot = 0
c   lpa = T means fit for disk pa
c   leps = T means fit for disk eps
c   lcenter = T means fit for disk xcen, ycen
      if( lpa )nd = nd + 1
      if( leps )nd = nd + 1
      if( lcentre )nd = nd + 2
      if( lvels )then
c assume that you always want to fit rotation velocity amplitudes
        ntot = nellip
c radial velocities add rows to linear matrix
c   force radial vel. ellipses to be same as circ. vel. ellipses
c   starting from nminr (.ge. 1) to nmaxr (.le. nellip)
        if( lradial )then
          nradial = nmaxr - nminr + 1
          ntot = ntot + nradial
        end if
c non-axisymmetric distortions add rows to linear matrix (nnasymm ellipses
c   for tangential part, nnasymm ellipses for radial part) as well as
c   non-linear terms
        if( lnax )then
c lphib = T means fit for bar angle
          if( lphib )then
            nd = nd + 1
          end if
          nnasymm = nmaxr - nminr + 1
          ntot = ntot + 2 * nnasymm
        end if
c systemic velocity adds a single row to linear matrix
        if( lsystemic )ntot = ntot + 1
c warp has 3 parameters
        if( lwarp )then
          if( lrwarp )nd = nd + 1
          if( lwepsm )nd = nd + 1
          if( lwpm )nd = nd + 1
          if( lradial .or. lnax )then
            print *, 'Both warp and non-axisymmetric flow selected'
            call crash( 'initlz' )
          end if
        end if
      else if( lphot )then
        if( ldisk )ntot = nellip
c fit a bar component - pa and ellipticity
        if( lnax )then
          if( lphib )nd = nd + 1
          if( lepsb )nd = nd + 1
          nbar = nmaxr - nminr + 1
          ntot = ntot + nbar
        else
          nbar = 0
        end if
c fit a parametric Sersic bulge component - flattening, index and r_eff
        if( lbulge )then
          if( lbleps )nd = nd + 1
          if( lsersn )nd = nd + 1
          if( lr_e )nd = nd + 1
          nbulge = mbulge
          if( .not. lI_0 )nbulge = 0
          ntot = ntot + nbulge
        else
          nbulge = 0
        end if
      else
        print *, 'Neither vels nor phot selected!'
        call crash( 'initlz' )
      end if
c check space
      if( nd .gt. md )then
        print *, 'Parameter md too small - value needed =', nd
        print *, 'Contact code administrator'
        call crash( 'initlz' )
      end if
      if( ntot .gt. mtot )then
        print *, 'Parameter mtot too small - value needed =', ntot
        print *, 'Contact code administrator'
        call crash( 'initlz' )
      end if
      return
      end
