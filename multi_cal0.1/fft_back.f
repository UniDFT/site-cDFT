      SUBROUTINE fft_back(array,fftable,ffwork,nfft1,nfft2,nfft3
     &     ,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)

************************************************************************
*   Time-stamp: <99/10/27 16:13:24 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb 13 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  array(*),fftable(*),ffwork(*)
      INTEGER nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable
     &     ,nffwork

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER isign

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      isign = -1

      CALL pubz3d(isign,nfft1,nfft2,nfft3,array,nfftdim1,nfftdim2
     &     ,fftable,nfftable,ffwork,nffwork)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

