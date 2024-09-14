c******************************************************
c* THIS PROGRAM IS TO GENERATE THE INITIAL DENSITY   ** 
c* PROFILE OF OXYGEN AND HYROGEN IN INHOMOGENEOUS    **
c* SYSTEM.                                           **              
C******************************************************


      SUBROUTINE Compute_initial_rho
      IMPLICIT NONE

      include 'param.h'
      include 'system.h'
      include 'numbers.h'
      INCLUDE 'rho.h'

*-------------  LOCAL VARIABLES   --------------
      INTEGER i, j, k, n, m
      REAL*8 x_com, y_com, z_com, rr

c** test part
      integer nb_hist, nhist, icos, jcos, nphi, ipusi
      real*8 ihhx, ihhy, ihhz, cos1, phi1, pusi1

      real*8 OMx1, Omy1, Omz1, rx, ry, rz

      real*8 norm_r,delta_r, hrom
      parameter(nb_hist = 624, delta_r = 0.025d0)

      real*8,ALLOCATABLE:: hrO(:),hrH(:)

      character*40 density_oxygen, density_hydrogen
      integer error
c interpolation
      real*8 GG
      integer n1, n2, n3

      real*8 x1, y1, z1, x2, y2, z2
      real*8 norm_k, kx, ky, kz, k2

      real*8 nr, r_nm
      real*8 cos_theta, phi


cc************************ test code below ******************************8



cc input the exact density profile of oxygen and hydrogen around solute water

      density_oxygen = "../Ne_water/density/hr_oxygen.dat"
      density_hydrogen = "../Ne_water/density/hr_hydrogen.dat"

      write(*,*)
      write(*,*)" # INPUT EXACT SOLVENT DENSITY PROFILE."
      write(*,*) "  INPUT FILES: ",density_oxygen, density_hydrogen
      write(*,*) "  delta_r =",delta_r, "nb_hist =",nb_hist

      ALLOCATE(hrO(0:nb_hist), stat=error)
      IF(error.ne.0) STOP "ALLOCATION hrO FAIL !"

      ALLOCATE(hrH(0:nb_hist), stat=error)
      IF(error.ne.0) STOP "ALLOCATION hrH FAIL !"

      open(11,file=density_oxygen)
      open(12,file=density_hydrogen)

      DO nhist =1, nb_hist

           read(11,*) norm_r, hrom
           hrO(nhist) = hrom

           read(12,*) norm_r, hrom
           hrH(nhist) = hrom ! density of one hydrogen site 

       END DO

       close(11)
       close(12)



       DO  k=1,nfft3
        DO  j=1,nfft2
         DO  i=1,nfft1

           x_com = (i-1)*Lg/FLOAT(nfft1) ! solvent oxygen_x
           y_com = (j-1)*Lg/FLOAT(nfft2)
           z_com = (k-1)*Lg/FLOAT(nfft3)

           rx = x_com - x_mol(1)
           ry = y_com - y_mol(1)
           rz = z_com - z_mol(1)

           r_nm = dsqrt(rx**2 + ry**2 + rz**2) ! solvent-solute distance

**

            nhist = INT(r_nm/delta_r)+1   
            nr = r_nm/delta_r+1 - nhist

            IF(nhist.lt.nb_hist) THEN
             hrom = hrO(nhist)*(1.-nr) + hrO(nhist+1)*nr
            ELSE
             hrom = zero
            END IF

            rho_O_old(i,j,k) = hrom + 1.0d0


            IF(nhist.lt.561) THEN
             hrom = hrH(nhist)*(1.-nr) + hrH(nhist+1)*nr
            ELSE
             hrom = zero
            END IF

            rho_H_old(i,j,k) = hrom + 1.0

           END DO
          END DO
        END DO

        DEALLOCATE(hrO)
        DEALLOCATE(hrH)


       OPEN(45,file = "../Results/rho_O_Z_exact.dat")
       OPEN(46,file = "../Results/rho_H_Z_exact.dat")
       DO j = nf2 + 1,nfft2

       write(45,*)  (j-1-nf2)*lg/nfft2
     &, rho_O_old(nf1,nf2,j)

       write(46,*)  (j-1-nf2)*lg/nfft2
     &, rho_H_old(nf1,nf2,j)

       END DO
       CLOSE(45)
       CLOSE(46)



       write(*,*) "  # Generate initial density profile."
       write(*,*)

       RETURN
       END



