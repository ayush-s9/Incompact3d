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

module jet

   USE decomp_2d
   USE variables
   USE param
 
   IMPLICIT NONE
 
   integer :: FS
   character(len=100) :: fileformat
   character(len=1),parameter :: NL=char(10) !new line character

   LOGICAL :: initialising
 
   PRIVATE ! All functions/subroutines private by default
   PUBLIC :: init_jet, boundary_conditions_jet, postprocess_jet, &
             visu_jet, visu_jet_init
 
 contains

   subroutine boundary_conditions_jet (rho,ux,uy,uz,phi)
 
     USE param
     USE variables
     USE decomp_2d
 
     implicit none
 
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
     
     
     call inflow (rho,phi)
     if (initialising) then
         return
     endif

     call outflow (rho,ux,uy,uz,phi)
 
     return
   end subroutine boundary_conditions_jet
   !********************************************************************
   subroutine inflow (rho,phi)
 
     USE param
     USE variables
     USE decomp_2d
     USE ibm_param
     use dbg_schemes, only: sin_prec, tanh_prec, sqrt_prec 
     use var, only: ta1, tb1

     implicit none
 
     integer  :: i,j,k,is
     real(mytype) :: D, r, y, z
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
 
     D = one

      do j = 1, xsize(2)
         y = real(j + xstart(2) - 2,mytype)*dy - half * yly
         do k = 1, xsize(3)
            z = real(k + xstart(3) - 2,mytype)*dz - half * zlz
            r = sqrt_prec(y**2 + z**2) + 1.e-16_mytype
            
            ! Setting up mean inflow - Michalke Jet Profile is used by default
            bxx1(j,k) = (u1 - u2) * half &
                          * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D))) + u2
            bxy1(j,k) = zero
            bxz1(j,k) = zero

            if (ilmn) then
               if (.not.ilmn_solve_temp) then
                  rho(1, j, k, 1) = dens1 -  (dens1 - dens2) * half &
                        * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D))) 
               else
                  rho(1, j, k, 1) = dens1 -  (dens1 - dens2) * half &
                        * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D))) 
                  !setting the temperature according to the desired density distribution
                  ta1(1, j, k) = ((pressure0 / rho(1, j, k, 1)) - one) / ((dens1 / dens2) - one)
               end if
            else 
               rho(1, j, k, 1) = one
            end if

            if (iscalar/=0) then
               do is = 1, numscalar
                  if (.not.massfrac(is)) then
                     phi(1, j, k, is) = half * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D)))
                  else if (is/=primary_species) then
                     phi(1, j, k, is) = one - half * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D)))
                  endif
               enddo

               if (primary_species.gt.0) then
                  phi(1, j, k, primary_species) = one

                  do is = 1, numscalar
                     if (massfrac(is).and.(is.ne.primary_species)) then
                        phi(1, j, k, primary_species) = phi(1, j, k, primary_species) - phi(1, j, k, is)
                     endif
                  enddo
               endif
            endif
            
         enddo
      enddo
      
      if (ilmn.and.ilmn_solve_temp) then
         CALL calc_rho_eos(rho(1,:,:,1), ta1(1,:,:), phi(1,:,:,:), tb1(1,:,:), 1, xsize(2), xsize(3))
      endif

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
   subroutine outflow (rho,ux,uy,uz,phi)
 
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
        if (primary_species.gt.0) then
          phi(nx,:,:,primary_species) = one
          do is = 1, numscalar
             if (massfrac(is).and.(is.ne.primary_species)) then
                phi(nx,:,:,primary_species) = phi(nx,:,:,primary_species) - phi(nx,:,:,is)
             endif
          enddo
        endif
     endif

     if (ilmn) then
       if (.not.ilmn_solve_temp) then
          rho(nx, :, :, 1) = rho(nx, :, :, 1) &
               - cx * (rho(nx, :, :, 1) - rho(nx - 1, :, :, 1))
       else
          !! Compute temperature at j-1:j to form advection equation
          CALL calc_temp_eos(ta1(nx-1:nx,:,:),rho(nx-1:nx,:,:,1),phi(nx-1:nx,:,:,:),tb1(nx-1:nx,:,:),2,xsize(2),xsize(3))
 
          ta1(nx,:,:) = ta1(nx,:,:) - cx * (ta1(nx,:,:) - ta1(nx-1,:,:))
 
          !! Need to compute rho (on boundary)
          CALL calc_rho_eos(rho(nx,:,:,1), ta1(nx,:,:), phi(nx,:,:,:), tb1(nx,:,:),1, xsize(2), xsize(3))
       endif
    endif

     if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
        write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)
 
     return
   end subroutine outflow
   !********************************************************************
   subroutine init_jet (rho1,ux1,uy1,uz1,ep1,phi1)
 
     USE decomp_2d
     USE decomp_2d_io
     USE variables
     USE param
     USE MPI
     use dbg_schemes, only: exp_prec
 
     implicit none
 
     real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
     real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

     logical :: col_init
 
     real(mytype) :: y,um
     integer :: k,j,i,ii,is,code
 
     if (iscalar==1) then
 
        phi1(:,:,:,:) = zero !change as much as you want
 
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
     call boundary_conditions_jet (rho1,ux1,uy1,uz1,phi1)
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

                 rho1(i,j,k,1)=rho1(1,j,k,1)
                 phi1(i,j,k,:)=phi1(1,j,k,:)
              enddo
           enddo
        enddo
     endif
 
#ifdef DEBG
   if (nrank .eq. 0) write(*,*) '# init end ok'
#endif
 
     return
   end subroutine init_jet
   !********************************************************************
 
   !############################################################################
   subroutine postprocess_jet(ux1,uy1,uz1,phi1,ep1)
 
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
 
   end subroutine postprocess_jet
 
   subroutine visu_jet_init (visu_initialised)
 
     use decomp_2d, only : mytype
     use decomp_2d_io, only : decomp_2d_register_variable
     use visu, only : io_name, output2D
     
     implicit none
 
     logical, intent(out) :: visu_initialised
 
     call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
     call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)
 
     visu_initialised = .true.
     
   end subroutine visu_jet_init
   !############################################################################
   !!
   !!  SUBROUTINE: visu_jet
   !!      AUTHOR: FS
   !! DESCRIPTION: Performs jetinder-specific visualization
   !!
   !############################################################################
   subroutine visu_jet(ux1, uy1, uz1, pp3, phi1, ep1, num)
 
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
     character(len=32), intent(in) :: num
 
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
     call write_field(di1, ".", "vort", trim(num), flush = .true.) ! Reusing temporary array, force flush
 
     !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
     di1 = zero
     di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                   - td1(:,:,:)*tb1(:,:,:) &
                   - tg1(:,:,:)*tc1(:,:,:) &
                   - th1(:,:,:)*tf1(:,:,:)
     call write_field(di1, ".", "critq", trim(num), flush = .true.) ! Reusing temporary array, force flush
 
   end subroutine visu_jet
 
 end module jet
 