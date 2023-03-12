!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module user_sim

  USE decomp_2d
  USE variables
  USE param
  USE dbg_schemes
  use var, only : zero, half, one, two

  IMPLICIT NONE

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1,volSimps1
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  logical :: initializing

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_user, boundary_conditions_user, postprocess_user, visu_user, &
            visu_user_init, set_fluid_properties_user

contains

  subroutine init_user (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ierror,is,code
    integer, dimension (:), allocatable :: seed
    integer ::  isize,ii
    logical :: bdry_init

    if (iscalar==1) then
      phi1(:,:,:,:) = zero
    endif

    ux1=zero; uy1=zero; uz1=zero;

    ! Generation of initial noise

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
               ux1(i,j,k)=init_noise*(two*ux1(i,j,k)-one)
               uy1(i,j,k)=init_noise*(two*uy1(i,j,k)-one)
               uz1(i,j,k)=init_noise*(two*uz1(i,j,k)-one)
            enddo
         enddo
      enddo
    endif

    if ((iscalar==1).and.(iin.ne.0)) then 
      call random_number(phi1)
      do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              phi1(i,j,k,:)=init_noise*(two*phi1(i,j,k,:)-one)
           enddo
        enddo
     enddo
    endif

    initializing=.true.
    call inflow(phi1)
    initializing=.false.

    do j = 1, xsize(2)
      do k = 1, xsize(3)
         ux1(1, j, k) = bxx1(j, k)
         uy1(1, j, k) = bxy1(j, k)
         uz1(1, j, k) = bxz1(j, k)
      end do
    end do

    ! Whether to initialize the domain using the inflow boundary
    bdry_init=.true.
    if ( bdry_init ) then
      do k=1,xsize(3)
        do j=1,xsize(2)
           do i=2,xsize(1)
              ux1(i,j,k)=ux1(i,j,k) + ux1(1,j,k)
              uy1(i,j,k)=uy1(i,j,k) + uy1(1,j,k)
              uz1(i,j,k)=uz1(i,j,k) + uz1(1,j,k)
              if ( iscalar==1 ) then
                  phi1(i,j,k,:)=phi1(1,j,k,:)
              end if
           enddo
        enddo
      enddo
    end if


#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_user

  subroutine boundary_conditions_user (ux,uy,uz,phi,ep)

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    ! Assuming x-direction to be the primary flow direction

    IF (nclx1.EQ.2) THEN
      call inflow(phi)
    ENDIF
    IF (nclxn.EQ.2) THEN
      call outflow(ux,uy,uz,phi)
    ENDIF
    ! wall boundary conditions for y=0 and y=Ly
    IF (ncly1.EQ.2) THEN
       byx1(:,:) = zero
       byy1(:,:) = zero
       byz1(:,:) = zero
       ! zero normal gradient boundary condition for concentration at the wall 
       if ( (iscalar==1).and.(nclyS1==2).and.(xstart(2)==1) ) then
          phi(:,1,:,:)=phi(:,2,:,:)
       end if
    ENDIF
    IF (nclyn.EQ.2) THEN
       byxn(:,:) = zero
       byyn(:,:) = zero
       byzn(:,:) = zero
       ! zero normal gradient boundary condition for concentration at the wall 
       if ( (iscalar==1).and.(nclySn==2).and.(xend(2)==ny) ) then
          phi(:,xsize(2),:,:)=phi(:,xsize(2)-1,:,:)
       end if
    ENDIF

    ! IF (nclz1.EQ.2) THEN
    ! ENDIF
    ! IF (nclzn.EQ.2) THEN
    ! ENDIF

  end subroutine boundary_conditions_user

  subroutine inflow (phi)
     
    USE param
    USE variables
    USE decomp_2d
    USE ibm_param
    USE dbg_schemes, only: sin_prec, tanh_prec, sqrt_prec, exp_prec, log_prec
   ! use var, only: ta1, tb1
    USE quadpack_generic
    
    implicit none
  
    integer  :: j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(2),xsize(3)) :: phio
    !quantities required for evaluating integral for velocity distribution
    real(mytype) :: y,int_result,h,abserr,resabs,resasc
    real(mytype) :: y_0
    real(mytype), parameter :: epsabs = 1.0E-06
    real(mytype), parameter :: epsrel = 1.0E-06
    integer, parameter :: limit = 1000
    integer, parameter :: lenw = limit*4
    integer, parameter :: key = 6
    integer :: ier, iwork(limit), last, neval
    real(mytype) :: work(lenw)
    
    !real(mytype),external :: u_dist
     
    h = 0.6_mytype
    y_0 = real(1 + xstart(2) - 2,mytype)*dy

     do j = 1, xsize(2)
        y = real(j + xstart(2) - 2,mytype)*dy  
        !call qag(u_dist,y_0,y,epsabs,epsrel,key,int_result,abserr,neval,ier)
        !call qk15(u_dist,y_0,y,int_result,abserr,resabs,resasc)
        call dqag(u_dist,y_0,y,epsabs,epsrel,key,int_result,abserr,neval,ier,limit,lenw,last,iwork,work)
        ! if(iscalar==0) int_result = 4.0_mytype * (y - y*y)
        bxx1(j,:) = real(int_result,mytype)
        bxy1(j,:) = zero
        bxz1(j,:) = zero
        y = real(y,mytype)
        if ( iscalar == 1) then
           phi(1,j,:,1) = half * (one + tanh_prec(40._mytype * (y - h)))
        end if
     enddo

     ! Mean Inflow + Inflow Noise
     call random_number(bxo)
     call random_number(byo)
     call random_number(bzo)
     call random_number(phio)
     do k=1,xsize(3)
        do j=1,xsize(2)
           bxx1(j,k)=bxx1(j,k) + (two * bxo(j,k) - one) * inflow_noise
           bxy1(j,k)=bxy1(j,k) + (two * byo(j,k) - one) * inflow_noise
           bxz1(j,k)=bxz1(j,k) + (two * bzo(j,k) - one) * inflow_noise
           phi(1,j,k,:)=phi(1,j,k,:) + (two * phio(j,k) - one) * inflow_noise
        enddo
     enddo  

    contains

        real(mytype) function u_dist(s)

           use dbg_schemes, only: tanh_prec, exp_prec, log_prec

           implicit none

           real(mytype), intent(in) :: s
           real(mytype), parameter :: h_ = 0.6
           real(mytype), parameter :: m = 2.0
           real(mytype), parameter :: rey = 500
           real(mytype), parameter :: press_grad = -0.034396385204418449
           real(mytype), parameter :: coeff_1 = rey*press_grad
           real(mytype), parameter :: coeff_2 = 7.3039516543885812

           u_dist = (coeff_1 * s + coeff_2) & 
                       * exp_prec(-half * (one + tanh_prec(40.0_mytype * (s - h_)))*log_prec(m))
     
        end function u_dist 

  end subroutine inflow

    !********************************************************************
  subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
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

    if (iscalar==1) then
       if (u2==zero) then
          cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
       elseif (u2==one) then
          cx=uxmax1*gdt(itr)*udx
       elseif (u2==two) then
          cx=u2*gdt(itr)*udx    !works better
       else
          stop
       endif

       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
          enddo
       enddo
    endif

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
  end subroutine outflow
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
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

  end subroutine postprocess_user

  subroutine visu_user_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    !call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    !call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

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
   !  if (sync_vel_needed) then
   !    call transpose_x_to_y(ux1,ux2)
   !    call transpose_x_to_y(uy1,uy2)
   !    call transpose_x_to_y(uz1,uz2)
   !    call transpose_y_to_z(ux2,ux3)
   !    call transpose_y_to_z(uy2,uy3)
   !    call transpose_y_to_z(uz2,uz3)
   !    sync_vel_needed = .false.
   !  endif

   !  !x-derivatives
   !  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
   !  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
   !  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
   !  !y-derivatives
   !  call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
   !  call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
   !  call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
   !  !!z-derivatives
   !  call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
   !  call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
   !  call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
   !  !!all back to x-pencils
   !  call transpose_z_to_y(ta3,td2)
   !  call transpose_z_to_y(tb3,te2)
   !  call transpose_z_to_y(tc3,tf2)
   !  call transpose_y_to_x(td2,tg1)
   !  call transpose_y_to_x(te2,th1)
   !  call transpose_y_to_x(tf2,ti1)
   !  call transpose_y_to_x(ta2,td1)
   !  call transpose_y_to_x(tb2,te1)
   !  call transpose_y_to_x(tc2,tf1)
   !  !du/dx=ta1 du/dy=td1 and du/dz=tg1
   !  !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
   !  !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
   !  !VORTICITY FIELD
   !  di1 = zero
   !  di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
   !                  + (tg1(:,:,:)-tc1(:,:,:))**2 &
   !                  + (tb1(:,:,:)-td1(:,:,:))**2)
   !  !call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

   !  !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
   !  di1 = zero
   !  di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
   !                - td1(:,:,:)*tb1(:,:,:) &
   !                - tg1(:,:,:)*tc1(:,:,:) &
   !                - th1(:,:,:)*tf1(:,:,:)
   !  call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush

    if ( iscalar==1 ) then
      call set_fluid_properties_user(phi1,viscosity)
      call write_field(viscosity,".","mu",num)
   end if
  end subroutine visu_user

  subroutine set_fluid_properties_user(phi1,mu1)

    implicit none
    
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: mu1

    integer :: i,j,k
    real(mytype), parameter :: m = two

    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             mu1(i,j,k) = exp_prec(phi1(i,j,k,1) * log_prec(m))
          end do
       end do
    end do

  end subroutine set_fluid_properties_user

end module user_sim
