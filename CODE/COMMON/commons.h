c Copyright (C) 2015, Jerry Sellwood and Kristine Spekkens
c
c updated by RZS Mar 2009
c      read a second fits file with velocity errors
c      need common block to pass fits array sizes
c updated by JAS Oct 2009
c      rationalized some variable names
c      logical for data type - FITS or text file
c updated by JAS Mar 2011
c      logicals phot and vels to branch to separate code
c      parameters for photo fitting
c updated by JAS Jul 2011
c      logical and parameters for a simple warp model
c updated by KS Mar 2012
c      light fraction variables and bulge profile added
c      updated file name variables
c updated by JAS Jan 2015
c      converted to f90 so that the sizes of all large arrays
c      are declared dynamically to make the code more flexible
c      the parameters mapx, mapy, mpx, & mpy have been eliminated
c updated by JAS Sep 2015
c      redefined file variable names - inpfile, invfile > datfile,
c      inevfile -> errfile, and added mskfile
c
c the following statement makes the main allocatable arrays visible in
c      every routine where it appears - it must come first, even before
c      the implicit statement
      use alldata
c
      implicit none
c
c global parameters
c
c md is the maximum number of non-linear parameters
      integer md
      parameter ( md = 9 )
c maximum numbers of ellipses
      integer mbar, mbulge, mellip, mtot
      parameter ( mellip = 250, mbar = mellip, mbulge = 1 )
      parameter ( mtot = 5 * mellip )
c
c mk is the maximum number of weights for each pixel
      integer mk
      parameter ( mk = 60 )
c
c maximum 1D spread of blurring raster, in pixels, for seeing corrections
      integer mblur
      parameter ( mblur = 61 )
c
c intrinsic disk thickness - zero for now
      real q
      parameter ( q = 0. )
c
c common blocks
c
      logical l2D, lerrfile, lFITS, lmask, lphot, ltext, lvels, verbose
      common / data_type / lFITS, ltext, verbose, lphot, lvels, l2D,
     +                     lerrfile, lmask
c
c input and output file names:
      character*100 datfile, errfile, infile, mskfile, outroot
      character*100 outpfile, outmfile, outrfile
      common / galname / infile, datfile, errfile, outpfile,
     +                   outmfile, outrfile, outroot, mskfile
c
      integer nd, nellip, nradial, ntot
      integer nbar, nbulge, nnasymm, nmaxr, nminr, nunc, seed
      logical lbulge, lcentre, ldisk, leps, lepsb, linter0, lnax
      logical lpa, lphib, lradial, lseeing, lsystemic, luncert, lwarp
      logical lbleps, lI_0, lr_e, lsersn, lrwarp, lwepsm, lwpm
      common / model / nd, nellip, nradial, ntot, nnasymm,
     +          nmaxr, nminr, nbar, nbulge, seed, nunc, lpa, leps,
     +          lcentre, lsystemic, lradial, lnax, lphib, luncert,
     +          linter0, lbulge, ldisk, lseeing, lwarp, lepsb,
     +          lr_e, lsersn, lbleps, lI_0, lrwarp, lwepsm, lwpm
c
      real eps, pa, incl, phib, phibprime, minr, maxr, junc
      real bar_pa, bar_eps, bulge_el, bulge_n, r_bulge, bulge_I0
      common / epars / eps, pa, incl, phib, phibprime, minr, maxr,
     +     junc, bar_pa, bar_eps, bulge_el, bulge_n, r_bulge, bulge_I0
c
      integer order
      common / morder / order
c
      real neweps, newpa, newincl, newphib
      real eneweps, enewpa, enewincl, enewphib
      real newbar_pa, newbar_eps, enewbar_pa, enewbar_eps
      real newbulge_l, newbulge_n, newr_bulge
      real enewbulge_l, enewbulge_n, enewr_bulge
      real newrw, newpm, newem, enewrw, enewpm, enewem
      common / newepars / neweps, newpa, newincl, newphib,
     +                    eneweps, enewpa, enewincl, enewphib,
     +                   newbar_pa, newbar_eps, enewbar_pa, enewbar_eps,
     +                    newbulge_l, newbulge_n, newr_bulge,
     +                    enewbulge_l, enewbulge_n, enewr_bulge,
     +                    newrw, newpm, newem, enewrw, enewpm, enewem
c
      real barfrac, blgfrac, dskfrac, ebarfrac, eblgfrac, edskfrac
      common / light /  barfrac, blgfrac, dskfrac, 
     +                  ebarfrac, eblgfrac, edskfrac
c
      real xcen, ycen, vsys
      common / center / xcen, ycen, vsys
c
      real newxcen, newycen
      real enewxcen, enewycen
      common / newcenter / newxcen, newycen, enewxcen, enewycen
c
      real sma0, smaf, sma( mellip )
      common / nelips / sma0, smaf, sma
c
      integer lox, loy, xrange, yrange
      common / fimage / lox, loy, xrange, yrange
c
      integer nsizex, nsizey
      integer istepout
      real pixscale, regeps, regpa, regrad, reginc
      logical VELmps
      common / image / nsizex, nsizey, pixscale,
     +                 istepout, regrad, regpa, regeps, VELmps
      equivalence ( reginc, regeps )
c
      integer pix_ind, inp_pts, pixct
      common / npixs / pix_ind, inp_pts, pixct
c
      real rwarp, wepsm, wphim
      common / warpp / rwarp, wphim, wepsm
c
      real*8 wba( mellip ), wcp( mellip ), wel( mellip )
      real*8 wphi( mellip ), wsi( mellip ), wsp( mellip )
      common / warpd / wel, wphi, wba, wcp, wsp, wsi
c
      logical fronly
      real b_over_a, b_over_a2, bulgel, bulgen, cosphib, cos2phib
      real cphi, rbulge, sini, sinphib, sin2phib, sphi
      common / wgtset / b_over_a, sini, cphi, sphi, cosphib, sinphib,
     +                  cos2phib, sin2phib, b_over_a2,
     +                  bulgel, bulgen, rbulge, fronly
c
      real lambda1, lambda2, typcl
      common / penalty / lambda1, lambda2, typcl
c
      integer iblur, jblur, nblur
      integer ibl( mblur ), jbl( mblur )
      real wblur( mblur )
      common / seecor / iblur, jblur, nblur, ibl, jbl, wblur
c
      real fitval( mtot ), efitval( mtot ), id( mellip )
      real ibulge( mbulge ), blgprof( mellip ), eblgprof( mellip )
      real ivsys, eivsys, irad( mellip ), ibar( mbar )
      real ibirad( mellip ), ibitan( mellip )
      real eid( mellip ), eirad( mellip ), eibulge( mbulge )
      real eibitan( mellip ), eibirad( mellip ), eibar( mbar )
      common / ints / fitval, efitval, id, irad, ivsys, ibirad, ibitan,
     +                eid, eirad, eivsys, eibirad, eibitan, 
     +                blgprof, eblgprof
      equivalence ( ibar( 1 ), irad( 1 ) ),
     +            ( ibulge( 1 ), ibirad( 1 ) ),
     +            ( eibar( 1 ), eirad( 1 ) ),
     +            ( eibulge( 1 ), eibirad( 1 ) )
c
      real dof
      common / freed / dof
c
      integer iter
      common / piter / iter
c
      real ringpts( mtot )
      common / rings / ringpts
c
      real final_chi2
      integer final_iter
      common / done / final_chi2, final_iter
c
      real sky, skysig, gain
      common / skyval / sky, skysig, gain
c
      real eISM, errtol, dsee
      common / sols / eISM, dsee, errtol
c
      real*8 pi
      parameter ( pi = 3.141592653589793238462643d0 )
