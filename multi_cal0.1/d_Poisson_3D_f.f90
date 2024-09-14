!*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 2006-2011 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
!
!*******************************************************************************
!  Content:
!  Fortran-90 double precision example of solving 3D Poisson problem in a
!  parallelepiped domain using MKL Poisson Library
!
!*******************************************************************************

program Poisson_3D_double_precision

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson

implicit none

integer nx, ny, nz
parameter(nx=4, ny=4, nz=4)
double precision pi
parameter(pi=3.14159265358979324D0)

integer ix, iy, iz, i, stat
integer ipar(128)
double precision ax, bx, ay, by, az, bz, lx, ly, lz, hx, hy, hz, xi, yi, zi
double precision cx, cy, cz, c1
double precision dpar(13*(nx+ny)/2+9)
! Note that proper packing of data in right-hand side array f and boundary
! arrays bd_* is automatically provided by the following declarations of arrays
double precision f(nx+1,ny+1,nz+1), u(nx+1,ny+1,nz+1)
double precision bd_ax(ny+1,nz+1), bd_bx(ny+1,nz+1)
double precision bd_ay(nx+1,nz+1), bd_by(nx+1,nz+1)
double precision bd_az(nx+1,ny+1), bd_bz(nx+1,ny+1)
double precision q
type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle
character(6) BCtype

! Printing the header for the example
print *, ''
print *, ' Example of use of MKL Poisson Library'
print *, ' **********************************************'
print *, ''
print *, ' This example gives the solution of 3D Poisson problem'
print *, ' with the equation -u_xx-u_yy-u_zz=f(x,y,z), 0<x<1, 0<y<1, 0<z<1,'
print *, ' f(x,y,z)=(12*pi*pi)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z),'
print *, ' and with the following boundary conditions:'
print *, '  u(0,y,z)=u(1,y,z)=1 (Dirichlet boundary conditions),'
print *, ' -u_y(x,0,z)=-2.0*pi*sin(2*pi*x)*sin(2*pi*z) (Neumann boundary condition),'
print *, '  u_y(x,1,z)= 2.0*pi*sin(2*pi*x)*sin(2*pi*z) (Neumann boundary condition),'
print *, ' -u_z(x,y,0)=-2.0*pi*sin(2*pi*x)*sin(2*pi*y) (Neumann boundary condition),'
print *, '  u_z(x,y,1)= 2.0*pi*sin(2*pi*x)*sin(2*pi*y) (Neumann boundary condition).'
print *, ' -----------------------------------------------------------------------'
print *, ' In general, the error should be of order O(1.0/nx^2+1.0/ny^2+1.0/nz^2)'
print '(1x,a,I1)', ' For this example, the value of nx=ny=nz is ', nx
print *, ' The approximation error should be of order 0.5D+0, if everything is OK'
print *, ' -----------------------------------------------------------------------'
print *, '                      DOUBLE PRECISION COMPUTATIONS                     '
print *, ' ======================================================================='
print *, ''

! Defining the parallelepiped domain 0<x<1, 0<y<1, 0<z<1 for 3D Poisson Solver
ax=0.0D0
bx=1.0D0
ay=0.0D0
by=1.0D0
az=0.0D0
bz=1.0D0

!*******************************************************************************
! Setting the coefficient q to 0.
! Note that this is the way to use Helmholtz Solver to solve Poisson problem!
!*******************************************************************************
q=0.0D0

! Computing the mesh size hx in x-direction
lx=bx-ax
hx=lx/nx
! Computing the mesh size hy in y-direction
ly=by-ay
hy=ly/ny
! Computing the mesh size hx in z-direction
lz=bz-az
hz=lz/nz

! Filling in the values of the TRUE solution
! u(x,y,z)=sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)+1
! in the mesh points into the array u
! Filling in the right-hand side
! f(x,y,z)=(12*pi*pi)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
! in the mesh points into the array f.
! We choose the right-hand side to correspond to the TRUE solution of
! Poisson equation.
! Here we are using the mesh sizes hx, hy, and hz computed before to compute
! the coordinates (xi,yi,zi) of the mesh points
do iz=1,nz+1
   do iy=1,ny+1
      do ix=1,nx+1
         xi=hx*(ix-1)/lx
         yi=hy*(iy-1)/ly
         zi=hz*(iz-1)/lz

         cx=dsin(2*pi*xi)
         cy=dsin(2*pi*yi)
         cz=dsin(2*pi*zi)

         u(ix,iy,iz)=1.0D0*cx*cy*cz
         f(ix,iy,iz)=(12.0D0*pi**2)*u(ix,iy,iz)
         u(ix,iy,iz)=u(ix,iy,iz)+1.0D0
      enddo
   enddo
enddo

! Setting the type of the boundary conditions on each surface
! of the parallelepiped domain:
! On the boundary laying on the plane x=0(=ax) Dirichlet boundary condition
! will be used
! On the boundary laying on the plane x=1(=bx) Dirichlet boundary condition
! will be used
! On the boundary laying on the plane y=0(=ay) Neumann boundary condition will be used
! On the boundary laying on the plane y=1(=by) Neumann boundary condition will be used
! On the boundary laying on the plane z=0(=az) Neumann boundary condition will be used
! On the boundary laying on the plane z=1(=bz) Neumann boundary condition will be used
BCtype = 'DDNNNN'

! Setting the values of the boundary function G(x,y,z) that is equal to
! the TRUE solution in the mesh points laying on Dirichlet boundaries
do iy = 1,ny+1
   do iz = 1,nz+1
     bd_ax(iy,iz) = 1.0D0
     bd_bx(iy,iz) = 1.0D0
   enddo
enddo
! Setting the values of the boundary function g(x,y,z) that is equal to
! the normal derivative of the TRUE solution in the mesh points laying on
! Neumann boundaries
do ix = 1,nx+1
   do iz = 1,nz+1
      bd_ay(ix,iz) = -2.D0*pi*dsin(2*pi*(ix-1)/nx)*dsin(2*pi*(iz-1)/nz)
      bd_by(ix,iz) =  2.D0*pi*dsin(2*pi*(ix-1)/nx)*dsin(2*pi*(iz-1)/nz)
   enddo
enddo
do ix = 1,nx+1
   do iy = 1,ny+1
      bd_az(ix,iy) = -2.D0*pi*dsin(2*pi*(ix-1)/nx)*dsin(2*pi*(iy-1)/ny)
      bd_bz(ix,iy) =  2.D0*pi*dsin(2*pi*(ix-1)/nx)*dsin(2*pi*(iy-1)/ny)
   enddo
enddo

! Initializing ipar array to make it free from garbage
do i=1,128
   ipar(i)=0
enddo

! Initializing simple data structures of Poisson Library for 3D Poisson Solver
call d_init_Helmholtz_3D(ax, bx, ay, by, az, bz, nx, ny, nz, BCtype, q, ipar, dpar,&
                                                                                 stat)
if (stat.ne.0) goto 999

! Initializing complex data structures of Poisson Library for 3D Poisson Solver
! NOTE: Right-hand side f may be altered after the Commit step. If you want
! to keep it, you should save it in another memory location!
call d_commit_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle,&
                                                            yhandle, ipar, dpar, stat)
if (stat.ne.0) goto 999

! Computing the approximate solution of 3D Poisson problem
! NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz
! should not be changed
! between the Commit step and the subsequent call to the Solver routine!
! Otherwise the results may be wrong.
call d_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle,&
                                                            yhandle, ipar, dpar, stat)
if (stat.ne.0) goto 999

! Cleaning the memory used by xhandle and yhandle
call free_Helmholtz_3D(xhandle, yhandle, ipar, stat)
if (stat.ne.0) goto 999
! Now we can use xhandle and yhandle to solve another 3D Poisson problem

! Printing the results
write(*,9)  nx
write(*,10) ny
write(*,11) nz
! Watching the error in the plane x=hx
ix=2
c1 = 0.0
do iy=1,ny+1
   do iz=1,nz+1
      write(*,12) (ix-1)*hx, (iy-1)*hy, (iz-1)*hz, f(ix,iy,iz)-u(ix,iy,iz)
      if (dabs(f(ix,iy,iz)-u(ix,iy,iz)).ge.c1) c1 = dabs(f(ix,iy,iz)-u(ix,iy,iz))
   enddo
enddo
print *, ''

if (c1.ge.0.5D+0) then
   print *, 'The computed solution seems to be inaccurate.'
   goto 999
endif

! Free MKL memory if any was allocated
call mkl_freebuffers
! Success message to print if everything is OK
print *, ' Double precision 3D Poisson example has successfully PASSED'
print *, ' through all steps of computation!'
stop 0
! Failure message to print if something went wrong
999 print *, 'Double precision 3D Poisson example FAILED to compute the solution...'
    ! Free MKL memory if any was allocated
    call mkl_freebuffers
    stop 1

9     format(1x,'The number of mesh intervals in x-direction is nx=',I1)
10    format(1x,'The number of mesh intervals in y-direction is ny=',I1)
11    format(1x,'The number of mesh intervals in z-direction is nz=',I1)
12    format(1x,'In the mesh point (',F4.2,',',F4.2,',',F4.2,') the error between ',&
                'the computed and the true solution is equal to ', E10.3)

! End of the example code
end
