!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
!!!! Trial -- Parabolic velocity profile as inlet and initial condition with our own geomcomplex subroutine
module user_sim

   USE decomp_2d
   USE variables
   USE param
 
   IMPLICIT NONE
 
   integer :: FS
   character(len=100) :: fileformat
   character(len=1),parameter :: NL=char(10) !new line character

   LOGICAL :: initialising
 
   PRIVATE ! All functions/subroutines private by default
   PUBLIC :: init_user, boundary_conditions_user, postprocess_user, &
             visu_user, visu_user_init, set_fluid_properties_user, &
             geomcomplex_user
 
 contains

   ! subroutine geomcomplex_user(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,remp)

   !    use decomp_2d, only : mytype
   !    use param, only : one, two, ten
   !    use ibm_param
   !    use dbg_schemes, only: sqrt_prec

   !    implicit none

   !    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
   !    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
   !    real(mytype),dimension(ny) :: yp
   !    real(mytype)               :: dx
   !    real(mytype)               :: remp
   !    integer                    :: i,j,k
   !    real(mytype)               :: y,z,r,rads2,kcon
   !    real(mytype)               :: dist_axi 
      
   !    ! Initialize epsi
   !    epsi(:,:,:)=zero

   !    ! Non-dimensionalization done with respect to the radius.
   !    do j = 1, xsize(2)
   !       y = real(j + xstart(2) - 2,mytype)*dy - half * yly
   !       do k = 1, xsize(3)
   !           z = real(k + xstart(3) - 2,mytype)*dz - half * zlz
   !           r = sqrt_prec(y**2 + z**2)

   !           if ( r.ge.one ) then
   !              epsi(:,j,k) = remp
   !           end if
   !       enddo
   !    end do

   ! end subroutine geomcomplex_user

!  SUBROUTINE geomcomplex_user(epsi,nx,nxi,nxf,ny,nyi,nyf,nz,nzi,nzf,xp,yp,zp,remp)
   SUBROUTINE geomcomplex_user(epsi,nx,nxi,nxf,ny,nyi,nyf,nz,nzi,nzf,dx,yp,remp)
   !
   !********************************************************************
 
     use decomp_2d, only : mytype
     USE MPI
     use param, only : zero, one, two,yly,zlz
   !   use param, only : new_rec
     use ibm
 
     implicit none
 
     !integer                                         :: nxi,nxf,ny,nyi,nyf,nzi,nzf
     integer                                         :: nx,nxi,nxf,ny,nyi,nyf,nz,nzi,nzf
     real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
   !   real(mytype),dimension(nx)                      :: xp
     real(mytype)                                    :: dx
     real(mytype),dimension(ny)                      :: yp
   !   real(mytype),dimension(nz)                      :: zp
     real(mytype)                                    :: remp
     real(mytype)                                    :: r,ym,zm,yc,zc,ra,tol
     !LOCALS
     !real(mytype),dimension(ny,nz)                   :: eps2d
     real(mytype),allocatable,dimension(:,:)         :: eps2d
     integer                                         :: i,j,k,irank,code,jj
 
     !==== New Reconstruction (through transverse-yz directions periodicity)? ====
   !   new_rec=1   ! 0: standard | 1: new reconstruction
     !============================================================================
 
     epsi(:,:,:) = zero
     yc = yly / two
     zc = zlz / two
     tol=1e-15
     ra=one
     !====DEBUG
     !!$tol=1e-16 ! BUG! breaks geometry symmetry
 
     allocate(eps2d(ny,nz))
 
     if (nrank.eq.0) then
         !Epsilon matrix (2D)
         eps2d(:,:) = zero
         do k=1,nz
            !  zm=zp(k)-zc
             zm=real(k + xstart(3) - 2,mytype)*dz-zc
             do j=1,ny
               !   ym=yp(j)-yc
                 ym=real(j + xstart(2) - 2,mytype)*dy-yc
                 r=dsqrt(ym*ym+zm*zm)
               !   if (r.gt.ra.and.r.lt.rao) then
                 if (r.gt.ra) then
                    eps2d(j,k)=remp
                 elseif (abs(r-ra).lt.tol) then
                     eps2d(j,k)=remp
               !   elseif (abs(r-rao).lt.tol) then
               !       eps2d(j,k)=remp
                 endif
             enddo
         enddo
 
         !Correct singular points
         !y-direction
         do k=1,nz
             do j=2,ny-1
                 if (eps2d(j-1,k).eq.zero.and.&
                     eps2d(j  ,k).eq.remp.and.&
                     eps2d(j+1,k).eq.zero) then
                     !singular solid
                     eps2d(j,k)=zero
                 !!$!-----------------------------------------------    
                 !!$!Attention! Correction of singular fluid points
                 !!$!can generate anomalies in the eps raf matrix
                 !!$elseif (eps2d(j-1,k).eq.remp.and.&
                 !!$        eps2d(j  ,k).eq.zero.and.&
                 !!$        eps2d(j+1,k).eq.remp) then
                 !!$    !singular fluid
                 !!$    eps2d(j,k)=remp
                 !!$!-----------------------------------------------    
                 endif
             enddo
         enddo
         !z-direction
         do j=1,ny
             do k=2,nz-1
                 if (eps2d(j,k-1).eq.zero.and.&
                     eps2d(j,k  ).eq.remp.and.&
                     eps2d(j,k+1).eq.zero) then
                     !singular solid
                     eps2d(j,k)=zero
                 !!$!-----------------------------------------------    
                 !!$!Attention! Correction of singular fluid points
                 !!$!can generate anomalies in the eps raf matrix
                 !!$elseif (eps2d(j,k-1).eq.remp.and.&
                 !!$        eps2d(j,k  ).eq.zero.and.&
                 !!$        eps2d(j,k+1).eq.remp) then
                 !!$    !singular fluid
                 !!$    eps2d(j,k)=remp
                 !!$!-----------------------------------------------    
                 endif
             enddo
         enddo
         !!$!====DEBUG
         !!$do k=1,nz
         !!$    zm=zp(k)-zc
         !!$    do j=1,ny
         !!$        ym=yp(j)-yc
         !!$        r=dsqrt(ym*ym+zm*zm)
         !!$        if (eps2d(j,k).eq.remp) print*, '(eps2d)', ym,zm,r-rao
         !!$    enddo
         !!$enddo
         !!$call sleep(2)
         !!$!stop
         !!$!=========
     endif
     call MPI_BCAST(eps2d,ny*nz,real_type,0,MPI_COMM_WORLD,code)
 
     !Epsilon matrix (3D)
     do k=nzi,nzf
         do j=nyi,nyf
             do i=nxi,nxf
                 epsi(i,j,k)=eps2d(j,k)
             enddo
         enddo
     enddo
     !
     deallocate(eps2d)
     !
     !!$!V1 (OLD)==========================================================================
     !!$!==== New Reconstruction (through transverse-yz directions periodicity)? ====
     !!$new_rec=1   ! 0: standard | 1: new reconstruction
     !!$!============================================================================
     !!$epsi(:,:,:) = zero
     !!$yc = yly / two
     !!$zc = zlz / two
     !!$!====DEBUG
     !!$!!$tol=1e-16 ! bug!
     !!$tol=1e-15
     !!$!Epsilon matrix
     !!$do k=nzi,nzf
     !!$    !====DEBUG
     !!$    !zm=real(k-1,mytype)*dz-zc
     !!$    zm=zp(k)-zc
     !!$    do j=nyi,nyf
     !!$        !====DEBUG
     !!$        ym=yp(j)-yc
     !!$        !ym=real(j-1,mytype)*dy-yc
     !!$        r=dsqrt(ym*ym+zm*zm)
     !!$        do i=nxi,nxf
     !!$            !if (r.ge.ra.and.r.le.rao) then
     !!$            if (r.gt.ra.and.r.lt.rao) then
     !!$               epsi(i,j,k)=remp
     !!$            !elseif (r.eq.ra) then
     !!$            elseif (abs(r-ra).lt.tol) then
     !!$                epsi(i,j,k)=remp
     !!$            !elseif (r.eq.rao) then
     !!$            elseif (abs(r-rao).lt.tol) then
     !!$                epsi(i,j,k)=remp
     !!$            endif
     !!$        enddo
     !!$    enddo
     !!$enddo
     !!$!==================================================================================
     !
     return
   end SUBROUTINE geomcomplex_user


   subroutine boundary_conditions_user (ux,uy,uz,phi,ep)
 
     USE param
     USE variables
     USE decomp_2d
     use dbg_schemes, only: sqrt_prec, abs_prec
 
     implicit none
 
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
     real(mytype) :: dr,y,z,ctheta,stheta,r,denominator
     integer :: i,j,k,jm,km
     
     call inflow (phi)
     if (initialising) then
         return
     endif

     call outflow (ux,uy,uz,phi)

     ! We can specify the lateral boundary conditions here
     ! The scalar boundary condition used is:
     ! zero gradient at the pipe wall

     if (iscalar.eq.1) then
        dr = sqrt_prec(dy**2 + dz**2)
        do i = 1, xsize(1)
           do j = 1, xsize(2)
               y = real(j + xstart(2) - 2,mytype)*dy - half * yly
               do k = 1, xsize(3)
                   z = real(k + xstart(3) - 2,mytype)*dz - half * zlz
                   r = sqrt_prec(y**2 + z**2)

                   if ( (r.le.(one)).or.(r.ge.(one - dr)) ) then
                     
                     ctheta = z/r
                     stheta = y/r 
                     
                     if (y.ge.zero) then
                         if ( z.ge.0 ) then
                            phi(i,j,k,:) = (phi(i,j,k-1,:)*dy*ctheta) + (phi(i,j-1,k,:)*dz*stheta)
                            phi(i,j,k,:) = phi(i,j,k,:) / (dy*ctheta + dz*stheta)
                         else
                            phi(i,j,k,:) = (-phi(i,j,k+1,:)*dy*ctheta) + (phi(i,j-1,k,:)*dz*stheta)
                            phi(i,j,k,:) = phi(i,j,k,:) / (-dy*ctheta + dz*stheta)
                         end if
                     else
                        if ( z.ge.0 ) then
                           phi(i,j,k,:) = (phi(i,j,k-1,:)*dy*ctheta) + (-phi(i,j+1,k,:)*dz*stheta)
                           phi(i,j,k,:) = phi(i,j,k,:) / (dy*ctheta - dz*stheta)
                        else
                           phi(i,j,k,:) = (phi(i,j,k+1,:)*dy*ctheta) + (phi(i,j+1,k,:)*dz*stheta)
                           phi(i,j,k,:) = phi(i,j,k,:) / (dy*ctheta + dz*stheta)
                        end if                        
                     endif
                   end if
               enddo
           end do
        end do
     endif

 
     return
   end subroutine boundary_conditions_user
   !********************************************************************
   subroutine inflow (phi)
  
     USE param
     USE variables
     USE decomp_2d
     USE ibm_param
     use dbg_schemes, only: sin_prec, tanh_prec, sqrt_prec 
     !use var, only: ta1, tb1

     implicit none
 
     integer  :: i,j,k,is
     real(mytype) :: r, y, z, interp_vel, r_i, delta
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
     
      r_i = 0.7_mytype
      delta = 0.02_mytype

      do j = 1, xsize(2)
         y = real(j + xstart(2) - 2,mytype)*dy - half * yly
         do k = 1, xsize(3)
            z = real(k + xstart(3) - 2,mytype)*dz - half * zlz
            r = sqrt_prec(y**2 + z**2) + 1.e-16_mytype
            
            bxx1(j,k) = zero
            bxy1(j,k) = zero
            bxz1(j,k) = zero
            
            if ( iscalar==1 ) then
               phi(1,j,k,:) = half * (one + erf((r - r_i) / delta))
            end if
            if ( r.le.one ) then
                call interpolated_profile(r,interp_vel)
                bxx1(j,k) = interp_vel

               ! bxx1(j,k) = one - r**2
            end if
         enddo
      enddo

      ! Mean Inflow + Inflow Noise
      call random_number(bxo)
      call random_number(byo)
      call random_number(bzo)
      do k=1,xsize(3)
         do j=1,xsize(2)
            bxx1(j,k)=bxx1(j,k) + (two * bxo(j,k) - one) * inflow_noise
            bxy1(j,k)=bxy1(j,k) + (two * byo(j,k) - one) * inflow_noise
            bxz1(j,k)=bxz1(j,k) + (two * bzo(j,k) - one) * inflow_noise
         enddo
      enddo  

     return
   end subroutine inflow
   !********************************************************************
   subroutine outflow (ux,uy,uz,phi)
 
     USE param
     USE variables
     USE decomp_2d
     USE MPI
     USE ibm_param
     use var, only: ta1, tb1
 
     implicit none
 
     integer :: j,k,code,is
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
     real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1
 
     udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz
 
     uxmax=-1609._mytype
     uxmin=1609._mytype
     do k=1,xsize(3)
        do j=1,xsize(2)
           if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
           if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
        enddo
     enddo
 
     call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)
 
     if (u1 == zero) then
        cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
     elseif (u1 == one) then
        cx=uxmax1*gdt(itr)*udx
     elseif (u1 == two) then
        cx=u2*gdt(itr)*udx    !works better
     else
        cx=(half*(u1+u2))*gdt(itr)*udx
     endif
 
     do k=1,xsize(3)
        do j=1,xsize(2)
           bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
           bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
           bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
        enddo
     enddo
 
     if (iscalar.ne.0) then
        if (u2==zero) then
           cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
        elseif (u2==one) then
           cx=uxmax1*gdt(itr)*udx
        elseif (u2==two) then
           cx=u2*gdt(itr)*udx    !works better
        else
           stop
        endif
        
        phi(nx,:,:,:)=phi(nx,:,:,:)-cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
     endif

     if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
        write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)
 
     return
   end subroutine outflow
   !********************************************************************
   subroutine init_user (ux1,uy1,uz1,ep1,phi1)
 
     USE decomp_2d
     USE decomp_2d_io
     USE variables
     USE param
     USE MPI
     use dbg_schemes, only: exp_prec
 
     implicit none
 
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

     logical :: col_init
 
     real(mytype) :: y,um
     integer :: k,j,i,ii,is,code
 
     if (iscalar==1) then
 
        phi1(:,:,:,:) = one !change as much as you want
 
     endif
 
     ux1=zero; uy1=zero; uz1=zero
 
     if (iin.ne.0) then
        call system_clock(count=code)
        if (iin.eq.2) code=0
        call random_seed(size = ii)
        call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))
 
        call random_number(ux1)
        call random_number(uy1)
        call random_number(uz1)
 
        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=1,xsize(1)
                 ux1(i,j,k)=init_noise * (two * ux1(i,j,k) - one)
                 uy1(i,j,k)=init_noise * (two * uy1(i,j,k) - one)
                 uz1(i,j,k)=init_noise * (two * uz1(i,j,k) - one)
              enddo
           enddo
        enddo
     endif

     initialising = .true.
     call boundary_conditions_user(ux1,uy1,uz1,phi1,ep1)
     initialising = .false.
     
     do j = 1, xsize(2)
         do k = 1, xsize(3)
            ux1(1, j, k) = bxx1(j, k)
            uy1(1, j, k) = bxy1(j, k)
            uz1(1, j, k) = bxz1(j, k)
         end do
     end do

     !If the domain needs to be initialized as a column, set flag col_init = .true.
     col_init = .true.

     if (col_init) then 
        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=2,xsize(1)
                 ux1(i,j,k)=ux1(i,j,k) + ux1(1,j,k)
                 uy1(i,j,k)=uy1(i,j,k) + uy1(1,j,k)
                 uz1(i,j,k)=uz1(i,j,k) + uz1(1,j,k)

                 phi1(i,j,k,:)=phi1(1,j,k,:)
              enddo
           enddo
        enddo
     endif
 
#ifdef DEBG
   if (nrank .eq. 0) write(*,*) '# init end ok'
#endif
 
     return
   end subroutine init_user
   !********************************************************************
 
   !############################################################################
   subroutine postprocess_user(ux1,uy1,uz1,phi1,ep1)
 
     USE MPI
     USE decomp_2d
     USE decomp_2d_io
     USE var, only : uvisu
     USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
     USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
     USE ibm_param
     use dbg_schemes, only: sqrt_prec
     
     real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
     real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
 
   end subroutine postprocess_user
 
   subroutine visu_user_init (visu_initialised)
 
     use decomp_2d, only : mytype
     use decomp_2d_io, only : decomp_2d_register_variable
     use visu, only : io_name, output2D
     
     implicit none
 
     logical, intent(out) :: visu_initialised
 
     call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
     call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)
 
     visu_initialised = .true.
     
   end subroutine visu_user_init
   !############################################################################
   !!
   !!  SUBROUTINE: visu_user
   !!      AUTHOR: FS
   !! DESCRIPTION: Performs user-specific visualization
   !!
   !############################################################################
   subroutine visu_user(ux1, uy1, uz1, pp3, phi1, ep1, num)
 
     use var, only : ux2, uy2, uz2, ux3, uy3, uz3
     USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
     USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
     use var, ONLY : nxmsize, nymsize, nzmsize
     use visu, only : write_field
     use ibm_param, only : ubcx,ubcy,ubcz
 
     implicit none
 
     real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
     real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
     real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
     real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
     real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: viscosity

     integer, intent(in) :: num
 
     ! Write vorticity as an example of post processing
 
     ! Perform communications if needed
     if (sync_vel_needed) then
       call transpose_x_to_y(ux1,ux2)
       call transpose_x_to_y(uy1,uy2)
       call transpose_x_to_y(uz1,uz2)
       call transpose_y_to_z(ux2,ux3)
       call transpose_y_to_z(uy2,uy3)
       call transpose_y_to_z(uz2,uz3)
       sync_vel_needed = .false.
     endif
 
     !x-derivatives
     call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
     call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
     call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
     !y-derivatives
     call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
     call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
     call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
     !!z-derivatives
     call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
     call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
     call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
     !!all back to x-pencils
     call transpose_z_to_y(ta3,td2)
     call transpose_z_to_y(tb3,te2)
     call transpose_z_to_y(tc3,tf2)
     call transpose_y_to_x(td2,tg1)
     call transpose_y_to_x(te2,th1)
     call transpose_y_to_x(tf2,ti1)
     call transpose_y_to_x(ta2,td1)
     call transpose_y_to_x(tb2,te1)
     call transpose_y_to_x(tc2,tf1)
     !du/dx=ta1 du/dy=td1 and du/dz=tg1
     !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
     !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
     !VORTICITY FIELD
     di1 = zero
     di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                     + (tg1(:,:,:)-tc1(:,:,:))**2 &
                     + (tb1(:,:,:)-td1(:,:,:))**2)
   !   call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush
 
     !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
     di1 = zero
     di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                   - td1(:,:,:)*tb1(:,:,:) &
                   - tg1(:,:,:)*tc1(:,:,:) &
                   - th1(:,:,:)*tf1(:,:,:)
   !   call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush
    if ( iscalar==1 ) then
       call set_fluid_properties_user(phi1,viscosity)
       call write_field(viscosity,".","mu",num)
    end if

   end subroutine visu_user
 
   subroutine set_fluid_properties_user(phi1,mu1)

      USE dbg_schemes, only: sin_prec, tanh_prec, sqrt_prec, exp_prec, log_prec
      implicit none
      
      real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
      real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: mu1
  
      integer :: i,j,k
      real(mytype), parameter :: m = 45._mytype
  
      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               mu1(i,j,k) = exp_prec(phi1(i,j,k,1) * log_prec(m))
            end do
         end do
      end do
  
    end subroutine set_fluid_properties_user

    subroutine interpolated_profile(r,interp_vel)

      use mod_splines

      implicit none

      real(mytype), intent(in) :: r
      real(mytype), intent(out) :: interp_vel

      real(mytype), dimension(1001) :: u_i, r_i
      type(spline) :: sp

      u_i = [ 3.8677016354583253,&
      3.8677016354583253, 3.8676936789928540, 3.8676701014555590, 3.8676309084163987, 3.8675761001347730, &
      3.8675056767937712, 3.8674196384391677, 3.8673179850495099, 3.8672007166278983, 3.8670678331773432, &
      3.8669193347007629, 3.8667552212009864, 3.8665754926807505, 3.8663801491427012, 3.8661691905893942, &
      3.8659426170232933, 3.8657004284467726, 3.8654426248621143, 3.8651692062715104, 3.8648801726770610, &
      3.8645755240807764, 3.8642552604845748, 3.8639193818902853, 3.8635678882996438, 3.8632007797139347, &
      3.8628180561285093, 3.8624197175415040, 3.8620057639529457, 3.8615761953628600, 3.8611310117712723, &
      3.8606702131782091, 3.8601937995836955, 3.8597017709877566, 3.8591941273904187, 3.8586708687917057, &
      3.8581319951916440, 3.8575775065902582, 3.8570074029875725, 3.8564216843836121, 3.8558203507784015, &
      3.8552034021719654, 3.8545708385643276, 3.8539226599555128, 3.8532588663455449, 3.8525794577344481, &
      3.8518844341222458, 3.8511737955089624, 3.8504475418946207, 3.8497056732792445, 3.8489481896628575, &
      3.8481750910454826, 3.8473863774271431, 3.8465820488078615, 3.8457621051876614, 3.8449265465665650, &
      3.8440753729445949, 3.8432085843217743, 3.8423261806981244, 3.8414281620736683, 3.8405145284484283, &
      3.8395852798224257, 3.8386404161956826, 3.8376799375682209, 3.8367038439400623, 3.8357121353112285, &
      3.8347048116817399, 3.8336818730516189, 3.8326433194208862, 3.8315891507895627, 3.8305193671576694, &
      3.8294339685252270, 3.8283329548922560, 3.8272163262587773, 3.8260840826248113, 3.8249362239903779, &
      3.8237727503554977, 3.8225936617201901, 3.8213989580844760, 3.8201886394483742, 3.8189627058119049, &
      3.8177211571750873, 3.8164639935379410, 3.8151912149004854, 3.8139028212627393, 3.8125988126247221, &
      3.8112791889864530, 3.8099439503479502, 3.8085930967092327, 3.8072266280703189, 3.8058445444312272, &
      3.8044468457919760, 3.8030335321525834, 3.8016046035130677, 3.8001600598734466, 3.7986999012337379, &
      3.7972241275939593, 3.7957327389541287, 3.7942257353142632, 3.7927031166743799, 3.7911648830344968, &
      3.7896110343946301, 3.7880415707547974, 3.7864564921150152, 3.7848557984753000, 3.7832394898356689, &
      3.7816075661961381, 3.7799600275567240, 3.7782968739174425, 3.7766181052783101, 3.7749237216393428, &
      3.7732137230005560, 3.7714881093619659, 3.7697468807235874, 3.7679900370854367, 3.7662175784475291, &
      3.7644295048098790, 3.7626258161725028, 3.7608065125354142, 3.7589715938986066, 3.7571210602617859, &
      3.7552549116248595, 3.7533731479878272, 3.7514757693506913, 3.7495627757134531, 3.7476341670761131, &
      3.7456899434386735, 3.7437301048011351, 3.7417546511634994, 3.7397635825257671, 3.7377568988879402, &
      3.7357346002500194, 3.7336966866120060, 3.7316431579739016, 3.7295740143357068, 3.7274892556974231, &
      3.7253888820590513, 3.7232728934205932, 3.7211412897820493, 3.7189940711434213, 3.7168312375047097, &
      3.7146527888659158, 3.7124587252270409, 3.7102490465880864, 3.7080237529490523, 3.7057828443099408, &
      3.7035263206707518, 3.7012541820314873, 3.6989664283921480, 3.6966630597527348, 3.6943440761132487, &
      3.6920094774736909, 3.6896592638340620, 3.6872934351943631, 3.6849119915545954, 3.6825149329147591, &
      3.6801022592748560, 3.6776739706348867, 3.6752300669948519, 3.6727705483547526, 3.6702954147145896, &
      3.6678046660743640, 3.6652983024340764, 3.6627763237937274, 3.6602387301533188, 3.6576855215128501, &
      3.6551166978723231, 3.6525322592317377, 3.6499322055910954, 3.6473165369503970, 3.6446852533096430, &
      3.6420383546688337, 3.6393758410279706, 3.6366977123870536, 3.6340039687460841, 3.6312946101050625, &
      3.6285696364639897, 3.6258290478228656, 3.6230728441816917, 3.6203010255404684, 3.6175135918991961, &
      3.6147105432578757, 3.6118918796165076, 3.6090576009750928, 3.6062077073336312, 3.6033421986921237, &
      3.6004610750505708, 3.5975643364089729, 3.5946519827673309, 3.5917240141256452, 3.5887804304839164, &
      3.5858212318421443, 3.5828464182003303, 3.5798559895584745, 3.5768499459165768, 3.5738282872746385, &
      3.5707910136326593, 3.5677381249906404, 3.5646696213485818, 3.5615855027064836, 3.5584857690643470, &
      3.5553704204221712, 3.5522394567799576, 3.5490928781377056, 3.5459306844954166, 3.5427528758530902, &
      3.5395594522107268, 3.5363504135683268, 3.5331257599258907, 3.5298854912834186, 3.5266296076409103, &
      3.5233581089983663, 3.5200709953557876, 3.5167682667131732, 3.5134499230705245, 3.5101159644278401, &
      3.5067663907851223, 3.5034012021423693, 3.5000203984995828, 3.4966239798567615, 3.4932119462139077, &
      3.4897842975710218, 3.4863410339281042, 3.4828821552851545, 3.4794076616421732, 3.4759175529991611, &
      3.4724118293561177, 3.4688904907130436, 3.4653535370699391, 3.4618009684268043, 3.4582327847836392, &
      3.4546489861404446, 3.4510495724972206, 3.4474345438539671, 3.4438039002106846, 3.4401576415673736, &
      3.4364957679240336, 3.4328182792806654, 3.4291251756372692, 3.4254164569938448, 3.4216921233503927, &
      3.4179521747069135, 3.4141966110634070, 3.4104254324198733, 3.4066386387763128, 3.4028362301327260, &
      3.3990182064891128, 3.3951845678454733, 3.3913353142018079, 3.3874704455581170, 3.3835899619144003, &
      3.3796938632706586, 3.3757821496268914, 3.3718548209830996, 3.3679118773392829, 3.3639533186954416, &
      3.3599791450515761, 3.3559893564076866, 3.3519839527637734, 3.3479629341198360, 3.3439263004758755, &
      3.3398740518318912, 3.3358061881878838, 3.3317227095438535, 3.3276236158998005, 3.3235089072557247, &
      3.3193785836116270, 3.3152326449675060, 3.3110710913233636, 3.3068939226791993, 3.3027011390350132, &
      3.2984927403908051, 3.2942687267465760, 3.2900290981023255, 3.2857738544580539, 3.2815029958137614, &
      3.2772165221694478, 3.2729144335251141, 3.2685967298807599, 3.2642634112363851, 3.2599144775919902, &
      3.2555499289475751, 3.2511697653031408, 3.2467739866586864, 3.2423625930142124, 3.2379355843697186, &
      3.2334929607252061, 3.2290347220806743, 3.2245608684361233, 3.2200713997915535, 3.2155663161469654, &
      3.2110456175023585, 3.2065093038577328, 3.2019573752130888, 3.1973898315684268, 3.1928066729237470, &
      3.1882078992790488, 3.1835935106343323, 3.1789635069895987, 3.1743178883448473, 3.1696566547000784, &
      3.1649798060552925, 3.1602873424104887, 3.1555792637656679, 3.1508555701208305, 3.1461162614759752, &
      3.1413613378311038, 3.1365907991862150, 3.1318046455413104, 3.1270028768963885, 3.1221854932514508, &
      3.1173524946064961, 3.1125038809615257, 3.1076396523165393, 3.1027598086715362, 3.0978643500265175, &
      3.0929532763814827, 3.0880265877364317, 3.0830842840913655, 3.0781263654462832, 3.0731528318011860, &
      3.0681636831560732, 3.0631589195109452, 3.0581385408658024, 3.0531025472206439, 3.0480509385754710, &
      3.0429837149302834, 3.0379008762850805, 3.0328024226398633, 3.0276883539946313, 3.0225586703493850, &
      3.0174133717041243, 3.0122524580588492, 3.0070759294135598, 3.0018837857682561, 2.9966760271229385, &
      2.9914526534776069, 2.9862136648322610, 2.9809590611869017, 2.9756888425415284, 2.9704030088961417, &
      2.9651015602507411, 2.9597844966053271, 2.9544518179598995, 2.9491035243144585, 2.9437396156690046, &
      2.9383600920235371, 2.9329649533780566, 2.9275541997325631, 2.9221278310870562, 2.9166858474415367, &
      2.9112282487960046, 2.9057550351504595, 2.9002662065049019, 2.8947617628593316, 2.8892417042137484, &
      2.8837060305681530, 2.8781547419225451, 2.8725878382769250, 2.8670053196312923, 2.8614071859856476, &
      2.8557934373399907, 2.8501640736943221, 2.8445190950486410, 2.8388585014029482, 2.8331822927572436, &
      2.8274904691115270, 2.8217830304657987, 2.8160599768200592, 2.8103213081743075, 2.8045670245285446, &
      2.7987971258827704, 2.7930116122369841, 2.7872104835911871, 2.7813937399453788, 2.7755613812995588, &
      2.7697134076537280, 2.7638498190078860, 2.7579706153620331, 2.7520757967161691, 2.7461653630702942, &
      2.7402393144244086, 2.7342976507785122, 2.7283403721326045, 2.7223674784866869, 2.7163789698407586, &
      2.7103748461948189, 2.7043551075488694, 2.6983197539029091, 2.6922687852569389, 2.6862022016109575, &
      2.6801200029649666, 2.6740221893189648, 2.6679087606729537, 2.6617797170269322, 2.6556350583808999, &
      2.6494747847348572, 2.6432988960888055, 2.6371073924427435, 2.6309002737966725, 2.6246775401505906, &
      2.6184391915044993, 2.6121852278583986, 2.6059156492122870, 2.5996304555661669, 2.5933296469200364, &
      2.5870132232738969, 2.5806811846277475, 2.5743335309815891, 2.5679702623354212, 2.5615913786892435, &
      2.5551968800430567, 2.5487867663968609, 2.5423610377506556, 2.5359196941044408, 2.5294627354582166, &
      2.5229901618119843, 2.5165019731657425, 2.5099981695194926, 2.5034787508732315, 2.4969437172269626, &
      2.4903930685806848, 2.4838268049343983, 2.4772449262881029, 2.4706474326417984, 2.4640343239954858, &
      2.4574056003491642, 2.4507612617028340, 2.4441013080564953, 2.4374257394101480, 2.4307345557637920, &
      2.4240277571174276, 2.4173053434710545, 2.4105673148246733, 2.4038136711782840, 2.3970444125318862, &
      2.3902595388854797, 2.3834590502390651, 2.3766429465926424, 2.3698112279462116, 2.3629638942997726, &
      2.3561009456533255, 2.3492223820068703, 2.3423282033604074, 2.3354184097139359, 2.3284930010674567, &
      2.3215519774209694, 2.3145953387744744, 2.3076230851279718, 2.3006352164814610, 2.2936317328349425, &
      2.2866126341884163, 2.2795779205418825, 2.2725275918953409, 2.2654616482487917, 2.2583800896022348, &
      2.2512829159556702, 2.2441701273090979, 2.2370417236625189, 2.2298977050159317, 2.2227380713693372, &
      2.2155628227227351, 2.2083719590761262, 2.2011654804295091, 2.1939433867828853, 2.1867056781362542, &
      2.1794523544896158, 2.1721834158429698, 2.1648988621963170, 2.1575986935496569, 2.1502829099029896, &
      2.1429515112563156, 2.1356044976096342, 2.1282418689629456, 2.1208636253162503, 2.1134697666695481, &
      2.1060602930228387, 2.0986352043761225, 2.0911945007293991, 2.0837381820826693, 2.0762662484359327, &
      2.0687786997891893, 2.0612755361424386, 2.0537567574956817, 2.0462223638489174, 2.0386723552021473, &
      2.0311067315553699, 2.0235254929085862, 2.0159286392617957, 2.0083161706149988, 2.0006880869681951, &
      1.9930443883213855, 1.9853850746745683, 1.9777101460277453, 1.9700196023809158, 1.9623134437340797, &
      1.9545916700872374, 1.9468542814403891, 1.9391012777935337, 1.9313326591466722, 1.9235484254998054, &
      1.9157485768529310, 1.9079331132060509, 1.9001020345591655, 1.8922553409122724, 1.8843930322653739, &
      1.8765151086184693, 1.8686215699715583, 1.8607124163246409, 1.8527876476777179, 1.8448472640307874, &
      1.8368912653838525, 1.8289196517369115, 1.8209324230899642, 1.8129295794430114, 1.8049111207960524, &
      1.7968770471490869, 1.7888273585021159, 1.7807620548551393, 1.7726811362113508, 1.7645846025787528, &
      1.7564724539568564, 1.7483446903447251, 1.7402013117414530, 1.7320423181461519, 1.7238677095579655, &
      1.7156774859760553, 1.7074716473996139, 1.6992501938278521, 1.6910131252600127, 1.6827604416953552, &
      1.6744921431331719, 1.6662082295727718, 1.6579087010134965, 1.6495935574547047, 1.6412627988957873, &
      1.6329164253361523, 1.6245544367752398, 1.6161768332125075, 1.6077836146474449, 1.5993747810795591, &
      1.5909503325083891, 1.5825102689334909, 1.5740545903544534, 1.5655832967708818, 1.5570963881824145, &
      1.5485938645887063, 1.5400757259894446, 1.5315419723843338, 1.5229926037731105, 1.5144276201555287, &
      1.5058470215313746, 1.4972508079004516, 1.4886389792625945, 1.4800115356176571, 1.4713684769655235, &
      1.4627098033060957, 1.4540355146393089, 1.4453456109651139, 1.4366400922834945, 1.4279189585944516, &
      1.4191822098980187, 1.4104298461942457, 1.4016618674832153, 1.3928782737650272, 1.3840790650398134, &
      1.3752642413077227, 1.3664338025689371, 1.3575877488236545, 1.3487260800721059, 1.3398487963145393, &
      1.3309558975512348, 1.3220473837824898, 1.3131232550086336, 1.3041835112300131, 1.2952281524470068, &
      1.2862571786600121, 1.2772705898694565, 1.2682683860757855, 1.2592505672794774, 1.2502171334810259, &
      1.2411680846809596, 1.2321034208798225, 1.2230231420781910, 1.2139272482766592, 1.2048157394758516, &
      1.1956886156764139, 1.1865458768790185, 1.1773875230843616, 1.1682135542931640, 1.1590239705061720, &
      1.1498187717241541, 1.1405979579479077, 1.1313615291782526, 1.1221094854160323, 1.1128418266621161, &
      1.1035585529173986, 1.0942596641827973, 1.0849451604592568, 1.0756150417477452, 1.0662693080492542, &
      1.0569079593648025, 1.0475309956954311, 1.0381384170422083, 1.0287302234062228, 1.0193064147885946, &
      1.0098669911904619, 1.0004119526129920, 0.9909412990573744, 0.9814550305248246, 0.9719531470165823, &
      0.9624356485339113, 0.9529025350781014, 0.9433538066504665, 0.9337894632523466, 0.9242095048851027, &
      0.9146139315501229, 0.9050027432488220, 0.8953759399826359, 0.8857335908937427, 0.8760757962630492, &
      0.8664025081167649, 0.8567136801627707, 0.8470092711678484, 0.8372892449576800, 0.8275535704168475, &
      0.8178022214888343, 0.8080351771760234, 0.7982524215396986, 0.7884539437000443, 0.7786397378361448, &
      0.7688098031859856, 0.7589641440464522, 0.7491027697733305, 0.7392256947813070, 0.7293329385439683, &
      0.7194245255938023, 0.7095004855221964, 0.6995608529794390, 0.6896056676747188, 0.6796349743761247, &
      0.6696488229106469, 0.6596472681641750, 0.6496303700814994, 0.6395981936663109, 0.6295508089812015, &
      0.6194882911476629, 0.6094107203460868, 0.5993181818157660, 0.5892129869478281, 0.5790983517456247, &
      0.5689721726522731, 0.5588327634887205, 0.5486790164331584, 0.5385104020210412, 0.5283269691450672, &
      0.5181293450551971, 0.5079187353586363, 0.4976969240198523, 0.4874662733605571, 0.4772297240597251, &
      0.4669907951535748, 0.4567535840355874, 0.4465227664564885, 0.4363035965242658, 0.4261019067041516, &
      0.4159241078186405, 0.4057771890474713, 0.3956687179276456, 0.3856068403534091, 0.3756002805762702, &
      0.3656582912388173, 0.3557913056084585, 0.3460164194891449, 0.3363523493663703, 0.3268178130968452, &
      0.3174315299085143, 0.3082122204005393, 0.2991786065433157, 0.2903494116784564, 0.2817433605188078, &
      0.2733791791484339, 0.2652755950226317, 0.2574513369679168, 0.2499251351820372, 0.2427166946853554, &
      0.2358446868210140, 0.2293149794650157, 0.2231304377858525, 0.2172920933443341, 0.2117991440935992, &
      0.2066489543791048, 0.2018370549386363, 0.1973571429022981, 0.1932010817925222, 0.1893589015240603, &
      0.1858187984039912, 0.1825671351317139, 0.1795884407989538, 0.1768614510391492, 0.1743621969356078, &
      0.1720741755766246, 0.1699813571154656, 0.1680680759539469, 0.1663190307424389, 0.1647192843798624, &
      0.1632542640136923, 0.1619097610399541, 0.1606719311032269, 0.1595272940966410, 0.1584627989453933, &
      0.1574704718826741, 0.1565442891020986, 0.1556771309346696, 0.1548624045593763, 0.1540940440031942, &
      0.1533665101410851, 0.1526747906959968, 0.1520144002388633, 0.1513813801886050, 0.1507722988121283, &
      0.1501842512243260, 0.1496148483149086, 0.1490616248071139, 0.1485224180462548, 0.1479954532505866, &
      0.1474790992189286, 0.1469718683306643, 0.1464724165457415, 0.1459795434046723, 0.1454921920285329, &
      0.1450094491189640, 0.1445305449581702, 0.1440548534089206, 0.1435818919145484, 0.1431111646119702, &
      0.1426420923511583, 0.1421743762748063, 0.1417077445369140, 0.1412419459009965, 0.1407767497400841, &
      0.1403119460367225, 0.1398473453829725, 0.1393827789804106, 0.1389180986401281, 0.1384531767827319, &
      0.1379879064383440, 0.1375222012466020, 0.1370559954566585, 0.1365892439271814, 0.1361219221263541, &
      0.1356540261318751, 0.1351855726309583, 0.1347165989203327, 0.1342471623029437, 0.1337772066268989, &
      0.1333066543837914, 0.1328355024935873, 0.1323637480300238, 0.1318913882206098, 0.1314184204466249, &
      0.1309448422431206, 0.1304706512989194, 0.1299958454566151, 0.1295204227125729, 0.1290443812169291, &
      0.1285677192735916, 0.1280904353402394, 0.1276125280283229, 0.1271339961030637, 0.1266548384834547, &
      0.1261750542422603, 0.1256946426060160, 0.1252136029550286, 0.1247319348233763, 0.1242496378989086, &
      0.1237667120232462, 0.1232831571917812, 0.1227989735536769, 0.1223141614118681, 0.1218287212230607, &
      0.1213426535977320, 0.1208559593001305, 0.1203686392482760, 0.1198806945139599, 0.1193921263227446, &
      0.1189029360539638, 0.1184131252407225, 0.1179226955698973, 0.1174316476026869, 0.1169399761752857, &
      0.1164476804970527, 0.1159547605688620, 0.1154612163915856, 0.1149670479660923, 0.1144722552932487, &
      0.1139768383739187, 0.1134807972089638, 0.1129841317992424, 0.1124868421456108, 0.1119889282489225, &
      0.1114903901100285, 0.1109912277297770, 0.1104914411090140, 0.1099910302485823, 0.1094899951493228, &
      0.1089883358120731, 0.1084860522376690, 0.1079831444269429, 0.1074796123807253, 0.1069754560998435, &
      0.1064706755851227, 0.1059652708373852, 0.1054592418574509, 0.1049525886461368, 0.1044453112042578, &
      0.1039374095326257, 0.1034288836320501, 0.1029197335033376, 0.1024099591472928, 0.1018995605647170, &
      0.1013885377564094, 0.1008768907231664, 0.1003646194657821, 0.0998517239850474, 0.0993382042817513, &
      0.0988240603566796, 0.0983092922106160, 0.0977938998443414, 0.0972778832586339, 0.0967612424542694, &
      0.0962439774320209, 0.0957260881926590, 0.0952075747369515, 0.0946884370656639, 0.0941686751795588, &
      0.0936482890793963, 0.0931272787659341, 0.0926056442399271, 0.0920833855021276, 0.0915605025532854, &
      0.0910369953941477, 0.0905128640254590, 0.0899881084479614, 0.0894627286623942, 0.0889367246694941, &
      0.0884100964699955, 0.0878828440646299, 0.0873549674541262, 0.0868264666392111, 0.0862973416206082, &
      0.0857675923990387, 0.0852372189752214, 0.0847062213498722, 0.0841745995237047, 0.0836423534974295, &
      0.0831094832717551, 0.0825759888473871, 0.0820418702250285, 0.0815071274053798, 0.0809717603891389, &
      0.0804357691770010, 0.0798991537696589, 0.0793619141678028, 0.0788240503721199, 0.0782855623832954, &
      0.0777464502020114, 0.0772067138289479, 0.0766663532647816, 0.0761253685101874, 0.0755837595658371, &
      0.0750415264324000, 0.0744986691105430, 0.0739551876009301, 0.0734110819042229, 0.0728663520210804, &
      0.0723209979521589, 0.0717750196981124, 0.0712284172595918, 0.0706811906372459, 0.0701333398317206, &
      0.0695848648436914, 0.0690357656738564, 0.0684860423228521, 0.0679356947913121, 0.0673847230798679, &
      0.0668331271891496, 0.0662809071197854, 0.0657280628724017, 0.0651745944476231, 0.0646205018460725, &
      0.0640657850683711, 0.0635104441151383, 0.0629544789869917, 0.0623978896845471, 0.0618406762084187, &
      0.0612828385592188, 0.0607243767375580, 0.0601652907440450, 0.0596055805792870, 0.0590452462438893, &
      0.0584842877384553, 0.0579227050635868, 0.0573604982198839, 0.0567976672079449, 0.0562342120283661, &
      0.0556701326817423, 0.0551054291686666, 0.0545401014897300, 0.0539741496455221, 0.0534075736366306, &
      0.0528403734636413, 0.0522725491271384, 0.0517041006277043, 0.0511350279659198, 0.0505653311423635, &
      0.0499950101576128, 0.0494240650122428, 0.0488524957068273, 0.0482803022419380, 0.0477074846181451, &
      0.0471340428360168, 0.0465599768961197, 0.0459852867990187, 0.0454099725452766, 0.0448340341354549, &
      0.0442574715701129, 0.0436802848498085, 0.0431024739750977, 0.0425240389465346, 0.0419449797646718, &
      0.0413652964300599, 0.0407849889432479, 0.0402040573047829, 0.0396225015152106, 0.0390403215750742, &
      0.0384575174849160, 0.0378740892452759, 0.0372900368566925, 0.0367053603197021, 0.0361200596348399, &
      0.0355341348026388, 0.0349475858236302, 0.0343604126983436, 0.0337726154273069, 0.0331841940110460, &
      0.0325951484500854, 0.0320054787449475, 0.0314151848961531, 0.0308242669042212, 0.0302327247696690, &
      0.0296405584930121, 0.0290477680747641, 0.0284543535154371, 0.0278603148155412, 0.0272656519755848, &
      0.0266703649960748, 0.0260744538775158, 0.0254779186204113, 0.0248807592252623, 0.0242829756925689, &
      0.0236845680228286, 0.0230855362165377, 0.0224858802741905, 0.0218856001962796, 0.0212846959832958, &
      0.0206831676357282, 0.0200810151540641, 0.0194782385387891, 0.0188748377903870, 0.0182708129093396, &
      0.0176661638961275, 0.0170608907512290, 0.0164549934751209, 0.0158484720682782, 0.0152413265311740, &
      0.0146335568642799, 0.0140251630680656, 0.0134161451429990, 0.0128065030895462, 0.0121962369081718, &
      0.0115853465993604, 0.0109738321636058, 0.0103616936013738, 0.0097489309131291, 0.0091355440993347, &
      0.0085215331604523, 0.0079068980969423, 0.0072916389092637, 0.0066757555978742, 0.0060592481632299, &
      0.0054421166057857, 0.0048243609259952, 0.0042059811243104, 0.0035869772011820, 0.0029673491570594, &
      0.0023470969923907, 0.0017262207076224, 0.0011047203031997, 0.0004825957795665, &
      0.0000000000000000]

      r_i = [ 0.0000000000000000,  &
      0.0000100000000000, 0.0010100000000000, 0.0020100000000000, 0.0030100000000000, 0.0040100000000000, &
      0.0050100000000000, 0.0060100000000000, 0.0070100000000000, 0.0080100000000000, 0.0090100000000000, &
      0.0100100000000000, 0.0110100000000000, 0.0120100000000000, 0.0130100000000000, 0.0140100000000000, &
      0.0150100000000000, 0.0160100000000000, 0.0170100000000000, 0.0180100000000000, 0.0190100000000000, &
      0.0200100000000000, 0.0210100000000000, 0.0220100000000000, 0.0230100000000000, 0.0240100000000000, &
      0.0250100000000000, 0.0260100000000000, 0.0270100000000000, 0.0280100000000000, 0.0290100000000000, &
      0.0300100000000000, 0.0310100000000000, 0.0320100000000000, 0.0330100000000000, 0.0340100000000000, &
      0.0350100000000000, 0.0360100000000000, 0.0370100000000000, 0.0380100000000000, 0.0390100000000000, &
      0.0400100000000000, 0.0410100000000000, 0.0420100000000000, 0.0430100000000000, 0.0440100000000000, &
      0.0450100000000000, 0.0460100000000000, 0.0470100000000000, 0.0480100000000000, 0.0490100000000000, &
      0.0500100000000000, 0.0510100000000000, 0.0520100000000000, 0.0530100000000000, 0.0540100000000000, &
      0.0550100000000000, 0.0560100000000000, 0.0570100000000000, 0.0580100000000000, 0.0590100000000000, &
      0.0600100000000000, 0.0610100000000000, 0.0620100000000000, 0.0630100000000000, 0.0640100000000000, &
      0.0650100000000000, 0.0660100000000000, 0.0670100000000000, 0.0680100000000000, 0.0690100000000000, &
      0.0700100000000000, 0.0710100000000000, 0.0720100000000000, 0.0730100000000000, 0.0740100000000000, &
      0.0750100000000000, 0.0760100000000000, 0.0770100000000000, 0.0780100000000000, 0.0790100000000000, &
      0.0800100000000000, 0.0810100000000000, 0.0820100000000000, 0.0830100000000000, 0.0840100000000000, &
      0.0850100000000000, 0.0860100000000000, 0.0870100000000000, 0.0880100000000000, 0.0890100000000000, &
      0.0900100000000000, 0.0910100000000000, 0.0920100000000000, 0.0930100000000000, 0.0940100000000000, &
      0.0950100000000000, 0.0960100000000000, 0.0970100000000000, 0.0980100000000000, 0.0990100000000000, &
      0.1000100000000000, 0.1010100000000000, 0.1020100000000000, 0.1030100000000000, 0.1040100000000000, &
      0.1050100000000000, 0.1060100000000000, 0.1070100000000000, 0.1080100000000000, 0.1090100000000000, &
      0.1100100000000000, 0.1110100000000000, 0.1120100000000000, 0.1130100000000000, 0.1140100000000000, &
      0.1150100000000000, 0.1160100000000000, 0.1170100000000000, 0.1180100000000000, 0.1190100000000000, &
      0.1200100000000000, 0.1210100000000000, 0.1220100000000000, 0.1230100000000000, 0.1240100000000000, &
      0.1250100000000000, 0.1260100000000000, 0.1270100000000000, 0.1280100000000000, 0.1290100000000000, &
      0.1300100000000000, 0.1310100000000000, 0.1320100000000000, 0.1330100000000000, 0.1340100000000000, &
      0.1350100000000000, 0.1360100000000000, 0.1370100000000000, 0.1380100000000000, 0.1390100000000000, &
      0.1400100000000000, 0.1410100000000000, 0.1420100000000000, 0.1430100000000000, 0.1440100000000000, &
      0.1450100000000000, 0.1460100000000000, 0.1470100000000000, 0.1480100000000000, 0.1490100000000000, &
      0.1500100000000000, 0.1510100000000000, 0.1520100000000000, 0.1530100000000000, 0.1540100000000000, &
      0.1550100000000000, 0.1560100000000000, 0.1570100000000000, 0.1580100000000000, 0.1590100000000000, &
      0.1600100000000000, 0.1610100000000000, 0.1620100000000000, 0.1630100000000000, 0.1640100000000000, &
      0.1650100000000000, 0.1660100000000000, 0.1670100000000000, 0.1680100000000000, 0.1690100000000000, &
      0.1700100000000000, 0.1710100000000000, 0.1720100000000000, 0.1730100000000000, 0.1740100000000000, &
      0.1750100000000000, 0.1760100000000000, 0.1770100000000000, 0.1780100000000000, 0.1790100000000000, &
      0.1800100000000000, 0.1810100000000000, 0.1820100000000000, 0.1830100000000000, 0.1840100000000000, &
      0.1850100000000000, 0.1860100000000000, 0.1870100000000000, 0.1880100000000000, 0.1890100000000000, &
      0.1900100000000000, 0.1910100000000000, 0.1920100000000000, 0.1930100000000000, 0.1940100000000000, &
      0.1950100000000000, 0.1960100000000000, 0.1970100000000000, 0.1980100000000000, 0.1990100000000000, &
      0.2000100000000000, 0.2010100000000000, 0.2020100000000000, 0.2030100000000000, 0.2040100000000000, &
      0.2050100000000000, 0.2060100000000000, 0.2070100000000000, 0.2080100000000000, 0.2090100000000000, &
      0.2100100000000000, 0.2110100000000000, 0.2120100000000000, 0.2130100000000000, 0.2140100000000000, &
      0.2150100000000000, 0.2160100000000000, 0.2170100000000000, 0.2180100000000000, 0.2190100000000000, &
      0.2200100000000000, 0.2210100000000000, 0.2220100000000000, 0.2230100000000000, 0.2240100000000000, &
      0.2250100000000000, 0.2260100000000000, 0.2270100000000000, 0.2280100000000000, 0.2290100000000000, &
      0.2300100000000000, 0.2310100000000000, 0.2320100000000000, 0.2330100000000000, 0.2340100000000000, &
      0.2350100000000000, 0.2360100000000000, 0.2370100000000000, 0.2380100000000000, 0.2390100000000000, &
      0.2400100000000000, 0.2410100000000000, 0.2420100000000000, 0.2430100000000000, 0.2440100000000000, &
      0.2450100000000000, 0.2460100000000000, 0.2470100000000000, 0.2480100000000000, 0.2490100000000000, &
      0.2500100000000000, 0.2510100000000000, 0.2520100000000000, 0.2530100000000000, 0.2540100000000000, &
      0.2550100000000000, 0.2560100000000000, 0.2570100000000000, 0.2580100000000000, 0.2590100000000000, &
      0.2600100000000000, 0.2610100000000000, 0.2620100000000000, 0.2630100000000000, 0.2640100000000000, &
      0.2650100000000000, 0.2660100000000000, 0.2670100000000000, 0.2680100000000000, 0.2690100000000000, &
      0.2700100000000000, 0.2710100000000000, 0.2720100000000000, 0.2730100000000000, 0.2740100000000000, &
      0.2750100000000000, 0.2760100000000000, 0.2770100000000000, 0.2780100000000000, 0.2790100000000000, &
      0.2800100000000000, 0.2810100000000000, 0.2820100000000000, 0.2830100000000000, 0.2840100000000000, &
      0.2850100000000000, 0.2860100000000000, 0.2870100000000000, 0.2880100000000000, 0.2890100000000000, &
      0.2900100000000000, 0.2910100000000000, 0.2920100000000000, 0.2930100000000000, 0.2940100000000000, &
      0.2950100000000000, 0.2960100000000000, 0.2970100000000000, 0.2980100000000000, 0.2990100000000000, &
      0.3000100000000000, 0.3010100000000000, 0.3020100000000000, 0.3030100000000000, 0.3040100000000000, &
      0.3050100000000000, 0.3060100000000000, 0.3070100000000000, 0.3080100000000000, 0.3090100000000000, &
      0.3100100000000000, 0.3110100000000000, 0.3120100000000000, 0.3130100000000000, 0.3140100000000000, &
      0.3150100000000000, 0.3160100000000000, 0.3170100000000000, 0.3180100000000000, 0.3190100000000000, &
      0.3200100000000000, 0.3210100000000000, 0.3220100000000000, 0.3230100000000000, 0.3240100000000000, &
      0.3250100000000000, 0.3260100000000000, 0.3270100000000000, 0.3280100000000000, 0.3290100000000000, &
      0.3300100000000000, 0.3310100000000000, 0.3320100000000000, 0.3330100000000000, 0.3340100000000000, &
      0.3350100000000000, 0.3360100000000000, 0.3370100000000000, 0.3380100000000000, 0.3390100000000000, &
      0.3400100000000000, 0.3410100000000000, 0.3420100000000000, 0.3430100000000000, 0.3440100000000000, &
      0.3450100000000000, 0.3460100000000000, 0.3470100000000000, 0.3480100000000000, 0.3490100000000000, &
      0.3500100000000000, 0.3510100000000000, 0.3520100000000000, 0.3530100000000000, 0.3540100000000000, &
      0.3550100000000000, 0.3560100000000000, 0.3570100000000000, 0.3580100000000000, 0.3590100000000000, &
      0.3600100000000000, 0.3610100000000000, 0.3620100000000000, 0.3630100000000000, 0.3640100000000000, &
      0.3650100000000000, 0.3660100000000000, 0.3670100000000000, 0.3680100000000000, 0.3690100000000000, &
      0.3700100000000000, 0.3710100000000000, 0.3720100000000000, 0.3730100000000000, 0.3740100000000000, &
      0.3750100000000000, 0.3760100000000000, 0.3770100000000000, 0.3780100000000000, 0.3790100000000000, &
      0.3800100000000000, 0.3810100000000000, 0.3820100000000000, 0.3830100000000000, 0.3840100000000000, &
      0.3850100000000000, 0.3860100000000000, 0.3870100000000000, 0.3880100000000000, 0.3890100000000000, &
      0.3900100000000000, 0.3910100000000000, 0.3920100000000000, 0.3930100000000000, 0.3940100000000000, &
      0.3950100000000000, 0.3960100000000000, 0.3970100000000000, 0.3980100000000000, 0.3990100000000000, &
      0.4000100000000000, 0.4010100000000000, 0.4020100000000000, 0.4030100000000000, 0.4040100000000000, &
      0.4050100000000000, 0.4060100000000000, 0.4070100000000000, 0.4080100000000000, 0.4090100000000000, &
      0.4100100000000000, 0.4110100000000000, 0.4120100000000000, 0.4130100000000000, 0.4140100000000000, &
      0.4150100000000000, 0.4160100000000000, 0.4170100000000000, 0.4180100000000000, 0.4190100000000000, &
      0.4200100000000000, 0.4210100000000000, 0.4220100000000000, 0.4230100000000000, 0.4240100000000000, &
      0.4250100000000000, 0.4260100000000000, 0.4270100000000000, 0.4280100000000000, 0.4290100000000000, &
      0.4300100000000000, 0.4310100000000000, 0.4320100000000000, 0.4330100000000000, 0.4340100000000000, &
      0.4350100000000000, 0.4360100000000000, 0.4370100000000000, 0.4380100000000000, 0.4390100000000000, &
      0.4400100000000000, 0.4410100000000000, 0.4420100000000000, 0.4430100000000000, 0.4440100000000000, &
      0.4450100000000000, 0.4460100000000000, 0.4470100000000000, 0.4480100000000000, 0.4490100000000000, &
      0.4500100000000000, 0.4510100000000000, 0.4520100000000000, 0.4530100000000000, 0.4540100000000000, &
      0.4550100000000000, 0.4560100000000000, 0.4570100000000000, 0.4580100000000000, 0.4590100000000000, &
      0.4600100000000000, 0.4610100000000000, 0.4620100000000000, 0.4630100000000000, 0.4640100000000000, &
      0.4650100000000000, 0.4660100000000000, 0.4670100000000000, 0.4680100000000000, 0.4690100000000000, &
      0.4700100000000000, 0.4710100000000000, 0.4720100000000000, 0.4730100000000000, 0.4740100000000000, &
      0.4750100000000000, 0.4760100000000000, 0.4770100000000000, 0.4780100000000000, 0.4790100000000000, &
      0.4800100000000000, 0.4810100000000000, 0.4820100000000000, 0.4830100000000000, 0.4840100000000000, &
      0.4850100000000000, 0.4860100000000000, 0.4870100000000000, 0.4880100000000000, 0.4890100000000000, &
      0.4900100000000000, 0.4910100000000000, 0.4920100000000000, 0.4930100000000000, 0.4940100000000000, &
      0.4950100000000000, 0.4960100000000000, 0.4970100000000000, 0.4980100000000000, 0.4990100000000000, &
      0.5000100000000000, 0.5010100000000000, 0.5020100000000000, 0.5030100000000000, 0.5040100000000000, &
      0.5050100000000000, 0.5060100000000000, 0.5070100000000000, 0.5080100000000000, 0.5090100000000000, &
      0.5100100000000000, 0.5110100000000000, 0.5120100000000000, 0.5130100000000000, 0.5140100000000000, &
      0.5150100000000000, 0.5160100000000000, 0.5170100000000000, 0.5180100000000000, 0.5190100000000000, &
      0.5200099999999999, 0.5210100000000000, 0.5220099999999999, 0.5230100000000000, 0.5240099999999999, &
      0.5250100000000000, 0.5260099999999999, 0.5270100000000000, 0.5280099999999999, 0.5290100000000000, &
      0.5300099999999999, 0.5310100000000000, 0.5320099999999999, 0.5330100000000000, 0.5340099999999999, &
      0.5350100000000000, 0.5360099999999999, 0.5370100000000000, 0.5380099999999999, 0.5390100000000000, &
      0.5400099999999999, 0.5410100000000000, 0.5420099999999999, 0.5430100000000000, 0.5440099999999999, &
      0.5450100000000000, 0.5460099999999999, 0.5470100000000000, 0.5480099999999999, 0.5490100000000000, &
      0.5500099999999999, 0.5510100000000000, 0.5520099999999999, 0.5530100000000000, 0.5540099999999999, &
      0.5550100000000000, 0.5560099999999999, 0.5570100000000000, 0.5580099999999999, 0.5590100000000000, &
      0.5600099999999999, 0.5610100000000000, 0.5620099999999999, 0.5630100000000000, 0.5640099999999999, &
      0.5650100000000000, 0.5660099999999999, 0.5670100000000000, 0.5680099999999999, 0.5690100000000000, &
      0.5700099999999999, 0.5710100000000000, 0.5720099999999999, 0.5730100000000000, 0.5740099999999999, &
      0.5750100000000000, 0.5760099999999999, 0.5770100000000000, 0.5780099999999999, 0.5790100000000000, &
      0.5800099999999999, 0.5810100000000000, 0.5820099999999999, 0.5830099999999999, 0.5840099999999999, &
      0.5850099999999999, 0.5860099999999999, 0.5870099999999999, 0.5880099999999999, 0.5890099999999999, &
      0.5900099999999999, 0.5910099999999999, 0.5920099999999999, 0.5930099999999999, 0.5940099999999999, &
      0.5950099999999999, 0.5960099999999999, 0.5970099999999999, 0.5980099999999999, 0.5990099999999999, &
      0.6000099999999999, 0.6010099999999999, 0.6020099999999999, 0.6030099999999999, 0.6040099999999999, &
      0.6050099999999999, 0.6060099999999999, 0.6070099999999999, 0.6080099999999999, 0.6090099999999999, &
      0.6100099999999999, 0.6110099999999999, 0.6120099999999999, 0.6130099999999999, 0.6140099999999999, &
      0.6150099999999999, 0.6160099999999999, 0.6170099999999999, 0.6180099999999999, 0.6190099999999999, &
      0.6200100000000000, 0.6210100000000000, 0.6220100000000000, 0.6230100000000000, 0.6240100000000000, &
      0.6250100000000000, 0.6260100000000000, 0.6270100000000000, 0.6280100000000000, 0.6290100000000000, &
      0.6300100000000000, 0.6310100000000000, 0.6320100000000000, 0.6330100000000000, 0.6340100000000000, &
      0.6350100000000000, 0.6360100000000000, 0.6370100000000000, 0.6380100000000000, 0.6390100000000000, &
      0.6400100000000000, 0.6410100000000000, 0.6420100000000000, 0.6430100000000000, 0.6440100000000000, &
      0.6450100000000000, 0.6460100000000000, 0.6470100000000000, 0.6480100000000000, 0.6490099999999999, &
      0.6500100000000000, 0.6510099999999999, 0.6520100000000000, 0.6530099999999999, 0.6540100000000000, &
      0.6550099999999999, 0.6560100000000000, 0.6570099999999999, 0.6580100000000000, 0.6590099999999999, &
      0.6600100000000000, 0.6610099999999999, 0.6620100000000000, 0.6630099999999999, 0.6640100000000000, &
      0.6650099999999999, 0.6660100000000000, 0.6670099999999999, 0.6680100000000000, 0.6690099999999999, &
      0.6700100000000000, 0.6710099999999999, 0.6720100000000000, 0.6730099999999999, 0.6740100000000000, &
      0.6750099999999999, 0.6760100000000000, 0.6770099999999999, 0.6780100000000000, 0.6790099999999999, &
      0.6800100000000000, 0.6810099999999999, 0.6820100000000000, 0.6830099999999999, 0.6840100000000000, &
      0.6850099999999999, 0.6860100000000000, 0.6870099999999999, 0.6880100000000000, 0.6890099999999999, &
      0.6900100000000000, 0.6910099999999999, 0.6920100000000000, 0.6930099999999999, 0.6940100000000000, &
      0.6950099999999999, 0.6960100000000000, 0.6970099999999999, 0.6980100000000000, 0.6990099999999999, &
      0.7000100000000000, 0.7010099999999999, 0.7020100000000000, 0.7030099999999999, 0.7040100000000000, &
      0.7050099999999999, 0.7060100000000000, 0.7070099999999999, 0.7080100000000000, 0.7090099999999999, &
      0.7100100000000000, 0.7110099999999999, 0.7120099999999999, 0.7130099999999999, 0.7140099999999999, &
      0.7150099999999999, 0.7160099999999999, 0.7170099999999999, 0.7180099999999999, 0.7190099999999999, &
      0.7200099999999999, 0.7210099999999999, 0.7220099999999999, 0.7230099999999999, 0.7240099999999999, &
      0.7250099999999999, 0.7260099999999999, 0.7270099999999999, 0.7280099999999999, 0.7290099999999999, &
      0.7300099999999999, 0.7310099999999999, 0.7320099999999999, 0.7330099999999999, 0.7340099999999999, &
      0.7350099999999999, 0.7360099999999999, 0.7370099999999999, 0.7380099999999999, 0.7390099999999999, &
      0.7400099999999999, 0.7410099999999999, 0.7420099999999999, 0.7430099999999999, 0.7440099999999999, &
      0.7450100000000000, 0.7460100000000000, 0.7470100000000000, 0.7480100000000000, 0.7490100000000000, &
      0.7500100000000000, 0.7510100000000000, 0.7520100000000000, 0.7530100000000000, 0.7540100000000000, &
      0.7550100000000000, 0.7560100000000000, 0.7570100000000000, 0.7580100000000000, 0.7590100000000000, &
      0.7600100000000000, 0.7610100000000000, 0.7620100000000000, 0.7630100000000000, 0.7640100000000000, &
      0.7650100000000000, 0.7660100000000000, 0.7670100000000000, 0.7680100000000000, 0.7690100000000000, &
      0.7700100000000000, 0.7710100000000000, 0.7720100000000000, 0.7730100000000000, 0.7740100000000000, &
      0.7750100000000000, 0.7760100000000000, 0.7770100000000000, 0.7780100000000000, 0.7790100000000000, &
      0.7800100000000000, 0.7810100000000000, 0.7820100000000000, 0.7830100000000000, 0.7840100000000000, &
      0.7850100000000000, 0.7860100000000000, 0.7870100000000000, 0.7880100000000000, 0.7890100000000000, &
      0.7900100000000000, 0.7910100000000000, 0.7920099999999999, 0.7930100000000000, 0.7940099999999999, &
      0.7950100000000000, 0.7960099999999999, 0.7970100000000000, 0.7980099999999999, 0.7990100000000000, &
      0.8000099999999999, 0.8010100000000000, 0.8020099999999999, 0.8030100000000000, 0.8040099999999999, &
      0.8050100000000000, 0.8060099999999999, 0.8070100000000000, 0.8080099999999999, 0.8090100000000000, &
      0.8100099999999999, 0.8110100000000000, 0.8120099999999999, 0.8130100000000000, 0.8140099999999999, &
      0.8150100000000000, 0.8160099999999999, 0.8170100000000000, 0.8180099999999999, 0.8190100000000000, &
      0.8200099999999999, 0.8210100000000000, 0.8220099999999999, 0.8230100000000000, 0.8240099999999999, &
      0.8250099999999999, 0.8260099999999999, 0.8270099999999999, 0.8280099999999999, 0.8290099999999999, &
      0.8300099999999999, 0.8310099999999999, 0.8320099999999999, 0.8330099999999999, 0.8340099999999999, &
      0.8350099999999999, 0.8360099999999999, 0.8370099999999999, 0.8380099999999999, 0.8390099999999999, &
      0.8400099999999999, 0.8410099999999999, 0.8420099999999999, 0.8430099999999999, 0.8440099999999999, &
      0.8450099999999999, 0.8460099999999999, 0.8470099999999999, 0.8480099999999999, 0.8490099999999999, &
      0.8500099999999999, 0.8510099999999999, 0.8520099999999999, 0.8530099999999999, 0.8540099999999999, &
      0.8550099999999999, 0.8560099999999999, 0.8570099999999999, 0.8580099999999999, 0.8590099999999999, &
      0.8600099999999999, 0.8610099999999999, 0.8620099999999999, 0.8630099999999999, 0.8640099999999999, &
      0.8650099999999999, 0.8660099999999999, 0.8670099999999999, 0.8680099999999999, 0.8690099999999999, &
      0.8700100000000000, 0.8710100000000000, 0.8720100000000000, 0.8730100000000000, 0.8740100000000000, &
      0.8750100000000000, 0.8760100000000000, 0.8770100000000000, 0.8780100000000000, 0.8790100000000000, &
      0.8800100000000000, 0.8810100000000000, 0.8820100000000000, 0.8830100000000000, 0.8840100000000000, &
      0.8850100000000000, 0.8860100000000000, 0.8870100000000000, 0.8880100000000000, 0.8890100000000000, &
      0.8900100000000000, 0.8910100000000000, 0.8920100000000000, 0.8930100000000000, 0.8940100000000000, &
      0.8950100000000000, 0.8960100000000000, 0.8970100000000000, 0.8980100000000000, 0.8990100000000000, &
      0.9000100000000000, 0.9010100000000000, 0.9020100000000000, 0.9030100000000000, 0.9040100000000000, &
      0.9050100000000000, 0.9060100000000000, 0.9070100000000000, 0.9080100000000000, 0.9090100000000000, &
      0.9100100000000000, 0.9110100000000000, 0.9120100000000000, 0.9130100000000000, 0.9140100000000000, &
      0.9150100000000000, 0.9160100000000000, 0.9170100000000000, 0.9180100000000000, 0.9190100000000000, &
      0.9200100000000000, 0.9210100000000000, 0.9220100000000000, 0.9230100000000000, 0.9240100000000000, &
      0.9250100000000000, 0.9260100000000000, 0.9270099999999999, 0.9280100000000000, 0.9290099999999999, &
      0.9300100000000000, 0.9310099999999999, 0.9320100000000000, 0.9330099999999999, 0.9340100000000000, &
      0.9350099999999999, 0.9360100000000000, 0.9370099999999999, 0.9380100000000000, 0.9390099999999999, &
      0.9400099999999999, 0.9410099999999999, 0.9420099999999999, 0.9430099999999999, 0.9440099999999999, &
      0.9450099999999999, 0.9460099999999999, 0.9470099999999999, 0.9480099999999999, 0.9490099999999999, &
      0.9500099999999999, 0.9510099999999999, 0.9520099999999999, 0.9530099999999999, 0.9540099999999999, &
      0.9550099999999999, 0.9560099999999999, 0.9570099999999999, 0.9580099999999999, 0.9590099999999999, &
      0.9600099999999999, 0.9610099999999999, 0.9620099999999999, 0.9630099999999999, 0.9640099999999999, &
      0.9650099999999999, 0.9660099999999999, 0.9670099999999999, 0.9680099999999999, 0.9690099999999999, &
      0.9700099999999999, 0.9710099999999999, 0.9720099999999999, 0.9730099999999999, 0.9740099999999999, &
      0.9750099999999999, 0.9760099999999999, 0.9770099999999999, 0.9780099999999999, 0.9790099999999999, &
      0.9800099999999999, 0.9810099999999999, 0.9820099999999999, 0.9830099999999999, 0.9840099999999999, &
      0.9850099999999999, 0.9860099999999999, 0.9870099999999999, 0.9880099999999999, 0.9890099999999999, &
      0.9900099999999999, 0.9910099999999999, 0.9920099999999999, 0.9930099999999999, 0.9940099999999999, &
      0.9950100000000000, 0.9960100000000000, 0.9970100000000000, 0.9980100000000000, &
      1.0000000000000000 ]

      sp = spline(r_i,u_i)
      interp_vel = sp%value(r)
      
    end subroutine interpolated_profile

 end module user_sim
 
