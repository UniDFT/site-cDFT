c************************************************************c
c* THIS SUBROUTINE IS TO CALCULATE THE EXTERNAL POTENTIAL  **
C* EXERTED ON EACH SITE DUE TO THE PRESENCE OF SOLUTE.     **
C*                              S.L. ZHAO                  **
C************************************************************


      SUBROUTINE Compute_Vext
      use vext
      implicit none

      include 'param.h'
      include 'fft.h'
      include 'system.h'
      include 'numbers.h'
      include 'Vext.h'

*-------------  LOCAL VARIABLES   --------------
      INTEGER i, j, k, n, m

      REAL*8 x_com, y_com, z_com, x_nm, y_nm, z_nm, x_m, y_m, z_m

      REAL*8 r_nm,r_nm2, r_nm6, sig_nm, sig_nm6, eps_nm, V_LJ, V_coul

      REAL*8 F,uljmin,uljtemp,blr

      real*8 rr2,rr

*BEGIN
      open(43, file ="../Results/VextO_Z.dat")
      open(44, file ="../Results/VextH_Z.dat")
c---- Compute the solute-solvent interaction ----

      DO  k=1,nfft3
        DO  j=1,nfft2
         DO  i=1,nfft1

** solvent site position
           x_com = (i-1)*Lg/FLOAT(nfft1)-lg/2.0
           y_com = (j-1)*Lg/FLOAT(nfft2)-lg/2.0
           z_com = (k-1)*Lg/FLOAT(nfft3)-lg/2.0
           rr2=x_com**2+y_com**2+z_com**2
           rr=sqrt(rr2)
           if(rr.gt.4.0)then
             VextO(i,j,k)=0.0D0
           else
             VextO(i,j,k)=5000.0D0
           endif
           VextH(i,j,k)=0.0D0
           VLJO(i,j,k) = VextO(i,j,k) 
          END DO 
          END DO 

cc  test code **
          write(43,*) (k-1)*lg/nfft3-lg/2.0, Beta*VextO(nf1,nf2,k)
          write(44,*) (k-1)*lg/nfft3-lg/2.0, Beta*VextH(nf1,nf2,k)
cc  test code **

          END DO  ! end the position on grid 
          close(43)
          close(44)

        write(*,*)"  # Calculate external potential."


       RETURN
       END





