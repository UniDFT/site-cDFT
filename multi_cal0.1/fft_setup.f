      SUBROUTINE fft_setup(fftable,nfft1,nfft2,nfft3,nfftable)

************************************************************************
*   Time-stamp: <99/10/27 16:16:35 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Tom Darden NIHS                                *
*              Modified by: Massimo MARCHI                             *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb 13 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  fftable(*)
      INTEGER nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      INTEGER nfftable

*------------------------- LOCAL VARIABLES ----------------------------*

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      write(6,*)'  # using public domain fft code'
      call pubz3di(nfft1,nfft2,nfft3,fftable,nfftable)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
