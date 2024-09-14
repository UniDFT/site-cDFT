      SUBROUTINE fft_backward(array,fftable,ffwork,nfft1,nfft2,nfft3
     &     ,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable
     &     ,nffwork
      REAL*8  array(2,nfftdim1,nfftdim2,nfftdim3)
      REAL*8  fftable(nfftable),ffwork(nffwork)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER isign
      INTEGER l,m,n

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      isign = -1

      CALL pubz3d(isign,nfft1,nfft2,nfft3,array,nfftdim1,nfftdim2
     &     ,fftable,nfftable,ffwork,nffwork)

*------------------------Normalization of back FFT --------------------*
      do l = 1, nfft1
       do m = 1,nfft2
         do n = 1,nfft3
            array(1,l,m,n) = array(1,l,m,n)/float(nfft1*nfft2*nfft3)
            array(2,l,m,n) = array(2,l,m,n)/float(nfft1*nfft2*nfft3)
         end do
       end do
      end do  
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

