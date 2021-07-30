      real function GHfit( x )
c  Copyright (C) 2019, Jerry Sellwood
c
c  function returns the value of the Gauss-Hermite fit to a line profile
      implicit none
c
c calling argument
      real x
c
c common blocks
c
      include 'comcube.h'
c
      real amp, h3, h4, mean, sigma, cont
      common / fitGH / amp, mean, sigma, h3, h4, cont
c
      real Hf3, Hf4, w
c
      w = ( x - mean ) / sigma
c
      Hf3 = w * ( 2 * w**2 - 3 ) / sqrt( 3. )
      Hf4 = ( 4 * w**4 - 12 * w**2 + 3 ) / sqrt( 24. )
c
      GHfit = cont + amp * exp( -0.5 * w**2 ) *
     +         ( 1 + h3 * Hf3 + h4 * Hf4 ) / ( sigma * sqrt( 2. * pi ) )
      return
      end
