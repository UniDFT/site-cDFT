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

      real*8 norm_r,delta_r, hrom, ratio
      parameter(nb_hist = 480, delta_r = 0.05d0)

      integer nb_cos1, nb_phi1, nb_pusi1
      parameter(nb_cos1 = 10, nb_phi1 = 20, nb_pusi1=10)
      REAL*8 delta_cos1, delta_phi1, delta_pusi1
      real*8,ALLOCATABLE:: hrO(:,:,:,:),hrH(:, :,:,:)

      character*40 density_oxygen, density_hydrogen
      integer error
c interpolation
      real*8 Deltacos(0:1), Deltaphi(0:1), Deltapusi(0:1)
      real*8 GG
      integer n1, n2, n3

      real*8 x1, y1, z1, x2, y2, z2
      real*8 norm_k, kx, ky, kz, k2

      real*8 nr, r_nm
      real*8 cos_theta, phi


cc************************ test code below ******************************8



cc input the exact density profile of oxygen and hydrogen around solute water
      delta_cos1 = two/nb_cos1
      delta_phi1 = twopi/nb_phi1
      delta_pusi1= pi/nb_pusi1

      density_oxygen = "../Water/hrom_relative_oxygen.dat"
      density_hydrogen = "../Water/hrom_relative_hydrogen.dat"

      write(*,*)
      write(*,*)" # INPUT EXACT SOLVENT DENSITY PROFILE."
      write(*,*) "  INPUT FILES: ",density_oxygen, density_hydrogen
      write(*,*) "  delta_r =",delta_r, "nb_hist =",nb_hist
      write(*,*) "  Orientations ",nb_cos1, nb_phi1, nb_pusi1

      ALLOCATE(hrO(0:nb_hist,nb_cos1,nb_phi1,nb_pusi1), stat=error)
      IF(error.ne.0) STOP "ALLOCATION hrO FAIL !"

      ALLOCATE(hrH(0:nb_hist,nb_cos1,nb_phi1,nb_pusi1), stat=error)
      IF(error.ne.0) STOP "ALLOCATION hrH FAIL !"

      open(11,file=density_oxygen)
      open(12,file=density_hydrogen)

!      ratio = 0.0332808332/0.033289065 
      ratio = 0.997

      DO nhist =0, nb_hist
        DO icos=1, nb_cos1
          DO nphi=1, nb_phi1
           DO ipusi =1, nb_pusi1

           read(11,*) norm_r, i, j, k, hrom
!           hrom = (hrom + 1.0d0)*float(4095)/float(4096) - 1.0d0
           hrom = (hrom + 1.0d0)*ratio - 1.0d0

           hrO(nhist,icos,nphi, ipusi) = hrom

           read(12,*) norm_r, i, j, k, hrom
!           hrom = (hrom + 1.0d0)*float(4095)/float(4096) - 1.0d0
           hrom = (hrom + 1.0d0)*ratio  - 1.0d0

           hrH(nhist,icos,nphi, ipusi) = hrom ! density of one hydrogen site 

          END DO
         END DO
        END DO

       END DO

       close(11)
       close(12)


      OMx1 = zero
      OMy1 = zero
      OMz1 = ONE    ! zero orientation of the solute molecule.

      ihhx = x_mol(2) - x_mol(3)
      ihhy = y_mol(2) - y_mol(3)
      ihhz = z_mol(2) - z_mol(3) ! H1H2 vector in solute

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

            

           CALL compute_angle(rx, ry, rz, cos_theta, phi) ! angle for solute-solvent vector

           call rotat(cos_theta,phi, OMx1,OMy1,OMz1, cos1,phi1)
c! solute orientation (cos1,phi1) within intermolecular frame


c-- we calculate pusi1 in the relative frame-------
          call transf_coord_fw(ihhx,ihhy,ihhz,cos_theta,phi,x1,y1,z1)
cc--calculate the projection of vector AB in the relative frame ^M

          

          call transf_coord_fw(x1,y1,z1,cos1,phi1,x2,y2,z2)
cc--calculate the projection of the vector in the molecule frame^M


          call compute_angle(x2,y2,z2,cos_theta,pusi1)
cc--the phi angular in the molcule frame is the pusi angular in the relative frame^M


cc!****  Now we are at the point: (r_nm, cos1, phi1, pusi1) *****

CC DENSITY PROFILE

            ipusi=int(pusi1/delta_pusi1) + 1
            IF(pusi1.lt.zero) STOP "pusi1 less than zero."
            if(ipusi.gt.nb_pusi1) ipusi=ipusi-nb_pusi1
            if(ipusi.gt.nb_pusi1) ipusi=ipusi-nb_pusi1


***
             icos = INT((one+cos1)/delta_cos1 + 0.5d0)
             if(icos.eq.0) icos=1
             if(icos.ge.nb_cos1) icos = nb_cos1 - 1

        Deltacos(1) = (cos1 + one)/delta_cos1 - (icos -0.5d0)
        Deltacos(0) = one - Deltacos(1)
*
               if(phi1.lt.zero) then
                  phi1= twopi + phi1
               end if
             nphi = INT(phi1/delta_phi1 + 0.5d0 )
             if(nphi.eq.0) nphi=1
             if(nphi.ge.nb_phi1) nphi = nb_phi1 - 1

        Deltaphi(1)  = phi1/delta_phi1 - (nphi-0.5d0)
        Deltaphi(0)  = one - Deltaphi(1)

**

            nhist = INT(r_nm/delta_r+0.5) !+1   
            nr = r_nm/delta_r - nhist

            IF(nhist.lt.nb_hist) THEN
!            IF(r_nm.lt.14.0) THEN  ! change

             hrom = zero
              do n1 = 0,1
               do n2 = 0,1

                   GG = hrO(nhist,icos+n1,nphi+n2,ipusi)
            hrom = hrom + GG * Deltacos(n1)*Deltaphi(n2)

                end do
               end do
                           
            ELSE
            hrom = zero
            END IF

            rho_O_old(i,j,k) = hrom + 1.0d0
            if(rho_O_old(i,j,k).lt.0.0D0)then
             rho_O_old(i,j,k)=0.0D0
            endif

            IF(nhist.lt.nb_hist) THEN
!            IF(r_nm.lt.14.0) THEN

             hrom = zero
              do n1 = 0,1
               do n2 = 0,1

                   GG = hrH(nhist,icos+n1,nphi+n2,ipusi)
            hrom = hrom + GG * Deltacos(n1)*Deltaphi(n2)

                end do
               end do

            ELSE
            hrom = zero
            END IF

            rho_H_old(i,j,k) = hrom + 1.0
            if(rho_H_old(i,j,k).lt.0.0D0)then
             rho_H_old(i,j,k)=0.0D0
            endif
           END DO
          END DO
        END DO

        DEALLOCATE(hrO)
        DEALLOCATE(hrH)


       OPEN(45,file = "../Results/rho_Or_exact.dat")
       OPEN(46,file = "../Results/rho_Hr_exact.dat")
       DO j = nf2 + 1,nfft2

       IF(mod(j,2).eq.1) then
       write(45,*)  (j-1-nf2)*lg/nfft2
     &, rho_O_old(nf1,j,nf3)
       write(46,*)  (j-1-nf2)*lg/nfft2
     &, rho_H_old(nf1,j,nf3)
       END IF

       END DO
       CLOSE(45)
       close(46)



       write(*,*) "  # Generate initial density profile."
       write(*,*)

       RETURN
       END



cc ###################################################################
cc this subroutine is to give the angular (theta,phi) of old vector (a,b,c)
cc in the new frame obtained by rotating the old frame theta1,phi1 angular 
cc ####################################################################
      subroutine rotat(cos_theta,phi,a,b,c,cos1,phi1)
      implicit none

      real*8 a,b,c,cos_theta,phi
      real*8 cos1,phi1
      real*8 xn,yn,zn

      call transf_coord_fw(a,b,c,cos_theta,phi,xn,yn,zn)

      call compute_angle(xn,yn,zn,cos1,phi1)

      return
       end


cc ####################################################################
cc this subroutine give the new coordinates (xn,yn,zn) in the new frame
cc givn we know the its coordinates (x,y,z) in the old frame
cc where new frame is obtained by rotating theta,phi of the old frame
cc ####################################################################
       subroutine transf_coord_fw(x,y,z,cos_theta,phi,xn,yn,zn)
       implicit none
       include 'numbers.h'

       real*8 x,y,z,cos_theta,phi
       real*8 xn,yn,zn
       real*8 sin_theta,sin_phi,cos_phi

       sin_theta=dsqrt(one-cos_theta**2)
       sin_phi=sin(phi)
       cos_phi=cos(phi)

      xn=x*cos_theta*cos_phi+y*cos_theta*sin_phi-z*sin_theta
      yn=-x*sin_phi+y*cos_phi
      zn=x*sin_theta*cos_phi+y*sin_theta*sin_phi+z*cos_theta

      return
      end

      subroutine compute_angle(x,y,z,cos_theta,phi)
      implicit none
      include 'numbers.h'

      real*8 x, y, z, cos_theta, phi
      real*8 rr,dd


      rr=sqrt(x**2+y**2+z**2)
      dd=sqrt(x**2+y**2)

cc-----------------------------
      if(dd.eq.zero) then
       phi=zero
      else if(y.ge.zero) then
       phi=acos(x/dd)
      else
       phi=twopi-acos(x/dd)
      end if

      if(phi.ge.twopi) phi=phi-twopi
cc-------- give the angular phi --------


      if(rr.eq.zero) then
       cos_theta=one
      else
       cos_theta= z/rr
      end if
cc ----------give the angular theta -----

       return
       end


      real*8 function  angle(x,y)
      implicit none
      include 'numbers.h'

      real*8 x,y,xx, r


        r=sqrt(x**2+y**2)

        if(r.eq.zero) then
            angle = zero
        else
           xx=x/r

         if(y.ge.0.0) then
             angle = acos(xx)
         else
             angle = twopi - acos(xx)
         end if

         end if

       return
       end

