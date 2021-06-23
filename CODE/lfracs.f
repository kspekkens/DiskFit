      subroutine lfracs( params )
c computes fractions of light in each fitted component
c
c   Created by JAS March 2012
c
      include 'commons.h'
c
c calling arguments
      real*8 params( md )
c
c local variables
      integer i, ij, ix, iy, j, k, km
      real barmod, bulgemod, diskmod, tmod
      real*8 barl, blgl, dskl, totl
c
c func will have skipped bad pixels unless seeing was turned on
c  need to reset weights in this case
      if( .not. lseeing )then
        call setwgt( params, .true., km )
      else
c need to redetermine the maximum number of ellipses contributing to any pixel
        km = 0
        do iy = 1, yrange
          do ix = 1, xrange
            do j = 1, ntot
              do i = 1, mk
                if( iwgt( i, ix, iy ) .eq. j )km = max( km, i )
              end do
            end do
          end do
        end do
      end if
c initialize in case any are never set
      diskmod = 0
      barmod = 0
      bulgemod = 0
      tmod = 0
c initialize total luminosities
      dskl = 0
      barl = 0
      blgl = 0
      totl = 0
c work over pixels
      do iy = 1, yrange
        do ix = 1, xrange
c select only pixels within the mask
          if( inmask( ix, iy ) )then
            tmod = 0
            k = 0
c sum model over components
            if( ldisk )then
              diskmod = 0
              do j = k + 1, k + nellip
                ij = 0
                do i = 1, km
                  if( iwgt( i, ix, iy ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )diskmod = diskmod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
              end do
              k = k + nellip
            end if
            tmod = diskmod
c bar
            if( lnax )then
              barmod = 0
              do j = k + 1, k + nbar
                ij = 0
                do i = 1, km
                  if( iwgt( i, ix, iy ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )barmod = barmod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
              end do
              k = k + nbar
              tmod = tmod + barmod
            end if
c bulge
            if( lbulge )then
              bulgemod = 0
              do j = k + 1, k + nbulge
                ij = 0
                do i = 1, km
                  if( iwgt( i, ix, iy ) .eq. j )ij = i
                end do
                if( ij .gt. 0 )bulgemod = bulgemod +
     +                                   wgt( ij, ix, iy ) * fitval( j )
              end do
              k = k + nbulge
              tmod = tmod + bulgemod
            end if
            totl = totl + tmod
            dskl = dskl + diskmod
            barl = barl + barmod
            blgl = blgl + bulgemod
          end if
        end do
      end do
c compute light fractions: note that diskfit always fits for a disk
      if( ldisk  )dskfrac = dskl / totl
      if( lnax   )barfrac = barl / totl
      if( lbulge )blgfrac = blgl / totl
      return
      end
