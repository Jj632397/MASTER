module mdl_explicit
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_rhs
  use mdl_result
  use mdl_bound
  use mdl_steps
  use mdl_time_delta
  use mdl_output
  use mdl_filter
  implicit none

  integer, private, save :: mode
  integer, private, save :: inner_iteration
  real(8), private, save :: relaxation_factor

contains

  subroutine explicit(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank
    integer :: isend_it,nstep_it
    integer :: nstep

    real(8), allocatable, dimension(:,:,:,:) :: QCR0
    real(8), allocatable, dimension(:,:,:,:) :: QCR1
    real(8), allocatable, dimension(:,:,:,:) :: QCR2
    real(8), allocatable, dimension(:,:,:,:) :: QCR3

    call mpisub_sbsp_get_myrank_world(myrank)

    call input_implicit(mode,inner_iteration,relaxation_factor)

    call steps_it_init

    call preconditioning_velocity_parameter(BLK)

    call time_delta(BLK)

    allocate(QCR0(BLK%imax_wh,BLK%jmax_wh,BLK%kmax,neqns))
    allocate(QCR1(BLK%imax_wh,BLK%jmax_wh,BLK%kmax,neqns))
    allocate(QCR2(BLK%imax_wh,BLK%jmax_wh,BLK%kmax,neqns))
    allocate(QCR3(BLK%imax_wh,BLK%jmax_wh,BLK%kmax,neqns))

    QCR0(:,:,:,:) = 0.0d0
    QCR1(:,:,:,:) = 0.0d0
    QCR2(:,:,:,:) = 0.0d0
    QCR3(:,:,:,:) = 0.0d0

    !Inner loop
    do while (1.eq.1) !*********************************************************

      QCR0(:,:,:,:) = BLK%QCS1(:,:,:,:)

      !1st-stage
      call right_hand_side(BLK)

      call explicit_func1

      BLK%QCS1(:,:,:,:) = QCR0(:,:,:,:) + 0.5d0*BLK%RHSS(:,:,:,:)

      if (i2dimension>0) then
        BLK%QCS1(:,:,:,4) = 0.0d0
      endif

      call result(BLK)

      call bound(BLK)

      call result(BLK)

      QCR1(:,:,:,:) = BLK%QCS1(:,:,:,:)

      !2nd-stage
      call right_hand_side(BLK)

      call explicit_func1

      BLK%QCS1(:,:,:,:) = QCR0(:,:,:,:) + 0.5d0*BLK%RHSS(:,:,:,:)

      if (i2dimension>0) then
        BLK%QCS1(:,:,:,4) = 0.0d0
      endif

      call result(BLK)

      call bound(BLK)

      call result(BLK)

      QCR2(:,:,:,:) = BLK%QCS1(:,:,:,:)

      !3rd-stage
      call right_hand_side(BLK)

      call explicit_func1

      BLK%QCS1(:,:,:,:) = QCR0(:,:,:,:) + BLK%RHSS(:,:,:,:)

      if (i2dimension>0) then
        BLK%QCS1(:,:,:,4) = 0.0d0
      endif

      call result(BLK)

      call bound(BLK)

      call result(BLK)

      QCR3(:,:,:,:) = BLK%QCS1(:,:,:,:)

      !4th-stage
      call right_hand_side(BLK)

      call explicit_func1

      BLK%QCS1(:,:,:,:) = QCR0(:,:,:,:) + QCR1(:,:,:,:)/6.0d0 + QCR2(:,:,:,:)/3.0d0 + QCR3(:,:,:,:)/3.0d0 + BLK%RHSS(:,:,:,:)/6.0d0

      if (i2dimension>0) then
        BLK%QCS1(:,:,:,4) = 0.0d0
      endif

      if (iact_lp_filter==1) then
        call filter(BLK)
      endif

      call result(BLK)

      call bound(BLK)

      call result(BLK)

      !****************************************************
      !Increment step count (Inner time iteration)        *
      call steps_it_increment !                           *
      !Output subroutine can be called hereafter.         *
      !****************************************************

      call steps_get_step_now(nstep)

      call steps_it_get_step_now(nstep_it)

      if (myrank==0) then
        open(34,file='data/progress.txt',form='formatted')
        write(34,*) nstep, nstep_it
        close(34)
      endif

      if (nstep==0.and.nstep_it==1) then
        call output_phys(BLK)
      endif

      if (nstep==0.and.nstep_it==5) then
        call output_phys(BLK)
      endif

      if (mod(nstep_it,5000)==0) then
        call output_phys(BLK)
      endif

      call output_l2norm_inner_iteration(BLK)

      call steps_it_isend(isend_it)
      if (isend_it.eq.1) exit
    enddo !*********************************************************************

    call explicit_func2

    call steps_it_final

    deallocate(QCR0)
    deallocate(QCR1)
    deallocate(QCR2)
    deallocate(QCR3)

  contains
!-------------------------------------------------------------------------------
    subroutine explicit_func1
      integer :: i,j,k
      real(8) :: cf1,cf2,cf3
      real(8) :: dt,rjb,bb,h,rh,u,v,w,hs,rrh,qq,c,c2
      real(8) :: hst,hsp,rht,rhp,hsqq,alph,dp,rdp,vp,vp2
      real(8) :: ys1,ys2,ys3,ys4,ys5,ys6,ys7
      real(8) :: rhb1,rhb2,rhb3,rhb4,rhb5,rhb6,rhb7,rhb8,rhb9,rhb10,rhb11,rhb12
      real(8) :: rhc1,rhc2,rhc3,rhc4,rhc5,rhc6,rhc7,rhc8,rhc9,rhc10,rhc11,rhc12

      cf1 = 1.5d0
      cf2 =-2.0d0
      cf3 = 0.5d0

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            dt  = BLK%DTLC(i,j,k)
            rjb = BLK%RJCB(i,j,k)

            rhb1 = BLK%RHSS(i,j,k,1)/rjb
            rhb2 = BLK%RHSS(i,j,k,2)/rjb
            rhb3 = BLK%RHSS(i,j,k,3)/rjb
            rhb4 = BLK%RHSS(i,j,k,4)/rjb
            rhb5 = BLK%RHSS(i,j,k,5)/rjb
            if (neqns>=6) rhb6 = BLK%RHSS(i,j,k, 6)/rjb
            if (neqns>=7) rhb7 = BLK%RHSS(i,j,k, 7)/rjb
            if (neqns>=8) rhb8 = BLK%RHSS(i,j,k, 8)/rjb
            if (neqns>=9) rhb9 = BLK%RHSS(i,j,k, 9)/rjb
            if (neqns>=10)rhb10= BLK%RHSS(i,j,k,10)/rjb
            if (neqns>=11)rhb11= BLK%RHSS(i,j,k,11)/rjb
            if (neqns>=12)rhb12= BLK%RHSS(i,j,k,12)/rjb

            bb = 1.0d0+1.5d0*dt/physical_time_step
            bb = 1.0d0/bb

            h = bb*dt

            rh= BLK%DNST(i,j,k)
            u = BLK%ULCT(i,j,k)
            v = BLK%VLCT(i,j,k)
            w = BLK%WLCT(i,j,k)
            hs= BLK%HTPS(i,j,k)
            if (neqns>=6) ys1 = BLK%YSPC(i,j,k,1)
            if (neqns>=7) ys2 = BLK%YSPC(i,j,k,2)
            if (neqns>=8) ys3 = BLK%YSPC(i,j,k,3)
            if (neqns>=9) ys4 = BLK%YSPC(i,j,k,4)
            if (neqns>=10)ys5 = BLK%YSPC(i,j,k,5)
            if (neqns>=11)ys6 = BLK%YSPC(i,j,k,6)
            if (neqns>=12)ys7 = BLK%YSPC(i,j,k,7)

            hst = BLK%SHCP(i,j,k)*CPREF
            hsp = BLK%DHPR(i,j,k)
            rht = BLK%DRTM(i,j,k)
            rhp = BLK%DRPR(i,j,k)

            rrh = 1.0d0/rh
            qq  = u*u+v*v+w*w

            hsqq = hs-0.5d0*qq

            alph = 1.0d0-rh*hsp

            dp = rh*rhp*hst+rht*alph
            rdp= 1.0d0/dp

            rhc1 = rdp*(alph-rhp*hsqq)*rhb1-rdp*rhp*u*rhb2-rdp*rhp*v*rhb3-rdp*rhp*w*rhb4+rdp*rhp*rhb5
            rhc2 =-rrh*u*rhb1+rrh*rhb2
            rhc3 =-rrh*v*rhb1+rrh*rhb3
            rhc4 =-rrh*w*rhb1+rrh*rhb4
            rhc5 = rdp*(rh*hst+rht*hsqq)*rhb1+rdp*rht*u*rhb2+rdp*rht*v*rhb3+rdp*rht*w*rhb4-rdp*rht*rhb5
            if (neqns>=6) rhc6 =-rrh*ys1*rhb1+rrh*rhb6
            if (neqns>=7) rhc7 =-rrh*ys2*rhb1+rrh*rhb7
            if (neqns>=8) rhc8 =-rrh*ys3*rhb1+rrh*rhb8
            if (neqns>=9) rhc9 =-rrh*ys4*rhb1+rrh*rhb9
            if (neqns>=10)rhc10=-rrh*ys5*rhb1+rrh*rhb10
            if (neqns>=11)rhc11=-rrh*ys6*rhb1+rrh*rhb11
            if (neqns>=12)rhc12=-rrh*ys7*rhb1+rrh*rhb12

            rhc1 = -(cf1*BLK%QCS1(i,j,k,1)+cf2*BLK%QCS2(i,j,k,1)+cf3*BLK%QCS3(i,j,k,1))/physical_time_step + rhc1
            rhc2 = -(cf1*BLK%QCS1(i,j,k,2)+cf2*BLK%QCS2(i,j,k,2)+cf3*BLK%QCS3(i,j,k,2))/physical_time_step + rhc2
            rhc3 = -(cf1*BLK%QCS1(i,j,k,3)+cf2*BLK%QCS2(i,j,k,3)+cf3*BLK%QCS3(i,j,k,3))/physical_time_step + rhc3
            rhc4 = -(cf1*BLK%QCS1(i,j,k,4)+cf2*BLK%QCS2(i,j,k,4)+cf3*BLK%QCS3(i,j,k,4))/physical_time_step + rhc4
            rhc5 = -(cf1*BLK%QCS1(i,j,k,5)+cf2*BLK%QCS2(i,j,k,5)+cf3*BLK%QCS3(i,j,k,5))/physical_time_step + rhc5
            if (neqns>=6) rhc6 = -(cf1*BLK%QCS1(i,j,k, 6)+cf2*BLK%QCS2(i,j,k, 6)+cf3*BLK%QCS3(i,j,k, 6))/physical_time_step + rhc6
            if (neqns>=7) rhc7 = -(cf1*BLK%QCS1(i,j,k, 7)+cf2*BLK%QCS2(i,j,k, 7)+cf3*BLK%QCS3(i,j,k, 7))/physical_time_step + rhc7
            if (neqns>=8) rhc8 = -(cf1*BLK%QCS1(i,j,k, 8)+cf2*BLK%QCS2(i,j,k, 8)+cf3*BLK%QCS3(i,j,k, 8))/physical_time_step + rhc8
            if (neqns>=9) rhc9 = -(cf1*BLK%QCS1(i,j,k, 9)+cf2*BLK%QCS2(i,j,k, 9)+cf3*BLK%QCS3(i,j,k, 9))/physical_time_step + rhc9
            if (neqns>=10)rhc10= -(cf1*BLK%QCS1(i,j,k,10)+cf2*BLK%QCS2(i,j,k,10)+cf3*BLK%QCS3(i,j,k,10))/physical_time_step + rhc10
            if (neqns>=11)rhc11= -(cf1*BLK%QCS1(i,j,k,11)+cf2*BLK%QCS2(i,j,k,11)+cf3*BLK%QCS3(i,j,k,11))/physical_time_step + rhc11
            if (neqns>=12)rhc12= -(cf1*BLK%QCS1(i,j,k,12)+cf2*BLK%QCS2(i,j,k,12)+cf3*BLK%QCS3(i,j,k,12))/physical_time_step + rhc12

            c  = BLK%CSOS(i,j,k)
            vp = dmin1(BLK%VPRE(i,j,k),c)

            c2 = c*c
            vp2= vp*vp

            rhb1 = rhc1-alph*(1.0d0-vp2/c2)/(rh*hst)*rhc5
            rhb2 = rhc2
            rhb3 = rhc3
            rhb4 = rhc4
            rhb5 = vp2/c2*rhc5
            if (neqns>=6) rhb6 = rhc6
            if (neqns>=7) rhb7 = rhc7
            if (neqns>=8) rhb8 = rhc8
            if (neqns>=9) rhb9 = rhc9
            if (neqns>=10)rhb10= rhc10
            if (neqns>=11)rhb11= rhc11
            if (neqns>=12)rhb12= rhc12

            BLK%RHSS(i,j,k,1) = h*rhb1
            BLK%RHSS(i,j,k,2) = h*rhb2
            BLK%RHSS(i,j,k,3) = h*rhb3
            BLK%RHSS(i,j,k,4) = h*rhb4
            BLK%RHSS(i,j,k,5) = h*rhb5
            if (neqns>=6) BLK%RHSS(i,j,k, 6) = h*rhb6
            if (neqns>=7) BLK%RHSS(i,j,k, 7) = h*rhb7
            if (neqns>=8) BLK%RHSS(i,j,k, 8) = h*rhb8
            if (neqns>=9) BLK%RHSS(i,j,k, 9) = h*rhb9
            if (neqns>=10)BLK%RHSS(i,j,k,10) = h*rhb10
            if (neqns>=11)BLK%RHSS(i,j,k,11) = h*rhb11
            if (neqns>=12)BLK%RHSS(i,j,k,12) = h*rhb12
          enddo
        enddo
      enddo

    end subroutine explicit_func1
!-------------------------------------------------------------------------------
    subroutine explicit_func2
      integer :: i,j,k,l
      real(8) :: qc1,qc2

      do l=1,neqns
        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            do i=BLK%ista,BLK%iend
              qc1 = BLK%QCS1(i,j,k,l)
              qc2 = BLK%QCS2(i,j,k,l)
              BLK%QCS2(i,j,k,l) = qc1
              BLK%QCS3(i,j,k,l) = qc2
            enddo
          enddo
        enddo
      enddo

    end subroutine explicit_func2

  end subroutine explicit
end module mdl_explicit
