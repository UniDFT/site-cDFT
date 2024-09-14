      subroutine get_fftdims(nfft1,nfft2,nfft3,nfftable)
      implicit none
      integer nfft1,nfft2,nfft3,nfftable
      integer n,nfftmax

      nfftmax = max(nfft1,nfft2,nfft3)
      nfftable = 4*nfftmax + 15

      return
      end
