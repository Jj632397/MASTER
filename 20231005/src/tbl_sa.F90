module mdl_tbl_sa
  use mdl_input_manager
  use mdl_mpisub_sbsp
  use mdl_decompo
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_fdm1
  use mdl_halo
  use mdl_steps
  implicit none

  real(8), public, parameter :: CB1SA = 0.1355d0
  real(8), public, parameter :: SGMSA = 2.0d0/3.0d0
  real(8), public, parameter :: CB2SA = 0.622d0
  real(8), public, parameter :: KPPSA = 0.41d0
  real(8), public, parameter :: CW2SA = 0.3d0
  real(8), public, parameter :: CW3SA = 2.0d0
  real(8), public, parameter :: CV1SA = 7.1d0
  real(8), public, parameter :: CT3SA = 1.2d0
  real(8), public, parameter :: CT4SA = 0.5d0
  real(8), public, parameter :: CW1SA = CB1SA/(KPPSA*KPPSA)+(1.0d0+CB2SA)/SGMSA

contains

  subroutine tbl_sa_result(BLK)
    type(block), intent(inout) :: BLK

    integer :: nstep,step_restart,mode_restart
    integer :: i,j,k
    real(8) :: vhat,rhd,visk,xai1,xai3,fv1,emt

    call steps_get_step_now(nstep)

    call input_restart(step_restart,mode_restart)

    if ((step_restart==0).and.(nstep<step_tbl_sa)) return

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          vhat = BLK%YSPC(i,j,k,neqn_tbl_sa-5)
          rhd  = BLK%DNST(i,j,k)*RHREF
          visk = BLK%VSSV(i,j,k)*VSREF/(BLK%DNST(i,j,k)*RHREF)

          xai1 = vhat/visk
          xai3 = xai1*xai1*xai1
          fv1 = xai3/(xai3+CV1SA*CV1SA*CV1SA)

          emt = rhd*vhat*fv1

          BLK%VSST(i,j,k) = emt/VSREF
          BLK%VSBT(i,j,k) = 0.0d0
          BLK%CKPT(i,j,k) = BLK%VSST(i,j,k)*BLK%SHCP(i,j,k)/tprntl
        enddo
      enddo
    enddo

  end subroutine tbl_sa_result
!-------------------------------------------------------------------------------
  subroutine tbl_sa_source(BLK,RHS)
    type(block), intent(inout) :: BLK
    real(8), intent(inout), dimension(:,:,:,:) :: RHS

    integer :: nstep,step_restart,mode_restart
    integer :: i,j,k
    real(8) :: rhx,rhy,rhz,uuy,uuz,vvx,vvz,wwx,wwy
    real(8) :: vsx,vsy,vsz
    real(8) :: vsvs,rhvs,vrtx,vhat,rhd,visk,dfw
    real(8) :: xai1,xai2,xai3,fv1,fv2,shat
    real(8) :: rsa,gsa,ft2,cw3sa6,fw
    real(8) :: cfdim,srct

    call steps_get_step_now(nstep)

    call input_restart(step_restart,mode_restart)

    if ((step_restart==0).and.(nstep<step_tbl_sa)) return

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          rhx = BLK%DDNX(i,j,k)
          rhy = BLK%DDNY(i,j,k)
          rhz = BLK%DDNZ(i,j,k)
          !--------------------
          uuy = BLK%DULY(i,j,k)
          uuz = BLK%DULZ(i,j,k)
          !--------------------
          vvx = BLK%DVLX(i,j,k)
          vvz = BLK%DVLZ(i,j,k)
          !--------------------
          wwx = BLK%DWLX(i,j,k)
          wwy = BLK%DWLY(i,j,k)
          !--------------------
          vsx = BLK%DYSX(i,j,k,neqn_tbl_sa-5)
          vsy = BLK%DYSY(i,j,k,neqn_tbl_sa-5)
          vsz = BLK%DYSZ(i,j,k,neqn_tbl_sa-5)

          vsvs = (vsx*vsx+vsy*vsy+vsz*vsz)/(XYREF*XYREF)
          rhvs = (rhx*vsx+rhy*vsy+rhz*vsz)*RHREF/(XYREF*XYREF)
          vrtx = dsqrt( (wwy-vvz)*(wwy-vvz) &
                       +(uuz-wwx)*(uuz-wwx) &
                       +(vvx-uuy)*(vvx-uuy) )*UUREF/XYREF

          vhat = BLK%YSPC(i,j,k,neqn_tbl_sa-5)
          rhd  = BLK%DNST(i,j,k)*RHREF
          visk = BLK%VSSV(i,j,k)*VSREF/(BLK%DNST(i,j,k)*RHREF)
          dfw  = BLK%YNKK(i,j,k)*XYREF

          dfw = dmax1(dfw,1.0d-16)

          xai1 = vhat/visk
          xai2 = xai1*xai1
          xai3 = xai2*xai1

          fv1 = xai3/(xai3+CV1SA*CV1SA*CV1SA)
          fv2 = 1.0d0-xai1/(1.0d0+xai1*fv1)
          shat= vrtx+vhat/(KPPSA*KPPSA*dfw*dfw)*fv2
          shat= DMAX1(shat,0.3d0*vrtx)

          rsa = DMIN1(vhat/(shat*KPPSA*KPPSA*dfw*dfw),10.0d0)
          gsa = rsa+CW2SA*(rsa**6-rsa)
          ft2 = CT3SA*DEXP(-CT4SA*xai2)
          cw3sa6 = CW3SA**6
          fw=gsa*((1.0d0+cw3sa6)/(gsa**6+cw3sa6))**(1.0d0/6.0d0)

          cfdim= XYREF/(RHREF*UUREF)
          srct = rhd*CB1SA*(1.0d0-ft2)*shat*vhat                      &
                -rhd*(CW1SA*fw-CB1SA*ft2/(KPPSA*KPPSA))*(vhat/dfw)**2 &
                +rhd*CB2SA*vsvs/SGMSA-(vhat+visk)*rhvs/SGMSA

          RHS(i,j,k,neqn_tbl_sa) = RHS(i,j,k,neqn_tbl_sa) + BLK%RJCB(i,j,k)*cfdim*srct
        enddo
      enddo
    enddo

  end subroutine tbl_sa_source

end module mdl_tbl_sa
