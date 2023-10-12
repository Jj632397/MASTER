module mdl_mean_and_rms
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_steps
  implicit none

  !input------------------------------------------------------------------------
  integer, private, parameter :: mode = 1 !<0>:inactive <1>:mean <2>:RMS(need mean)
  integer, private, parameter :: mode_Tpipe = 1 !<0>:inactive or usual mean <1>:mean at outlet !! FURUSAWA
  integer, private, parameter :: step_rst_mean = 100000
  !-----------------------------------------------------------------------------

  character(18), private, parameter :: prefix = "restart/rst_mean_n"
  character(4 ), private, parameter :: suffix = ".dat"

  integer, private, parameter :: lmax = 8

  real(8), save, private, allocatable, dimension(:,:,:,:) :: QMEN1
  real(8), save, private, allocatable, dimension(:,:,:,:) :: QMEN2
  real(8), save, private, allocatable, dimension(:,:,:,:) :: QRMS
  !1:Density
  !2:Velocity_u
  !3:Velocity_v
  !4:Velocity_w
  !5:Pressure
  !6:Temperature
  !7:Viscosity coefficient
  !8:Distance from wall

contains

  subroutine mean_and_rms(BLK)
    type(block), intent(inout) :: BLK

    if (mode==0) return

    if(mode_Tpipe==0) then  !! FURUSAWA
      
        
      if (mode==2) then
        call sub_rms(BLK)
      endif
    else if(mode_Tpipe==1) then  !! FURUSAWA
      call sub_mean_Tpipe(BLK)  !! FURUSAWA
      call sub_reaction_rate(BLK)  !! KENTA
      call sub_mean(BLK) !furukawa
    endif

  end subroutine mean_and_rms
!-------------------------------------------------------------------------------
  subroutine sub_mean_Tpipe(BLK)
  
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: nstep,step_start,step_final,nperiod
    integer :: i,j,k,ii,jj,kk
    real(8) :: rh,u,v,w,p,t,cmu,sc01,sc02
    real(8) :: rnperi
    character(100) :: fname

900 format('data/mean_n_t_pipe',i7.7,'.dat')

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call steps_get_step_now(nstep)
    call steps_get_step_start(step_start)
    call steps_get_step_final(step_final)

    nperiod = 25000

    if (nstep==step_start+1) then
      allocate(QMEN1(imax,jmax,kmax,lmax))
      QMEN1(:,:,:,:) = 0.0d0
    endif

    if (step_start+1<=nstep.and.nstep<=step_final) then

      rnperi = 1.0d0/dble(nperiod)

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1

            rh = BLK%DNST(i,j,k)  *RHREF
            u  = BLK%ULCT(i,j,k)  *UUREF
            v  = BLK%VLCT(i,j,k)  *UUREF
            w  = BLK%WLCT(i,j,k)  *UUREF
            p  = BLK%PRSS(i,j,k)  *PPREF
            t  = BLK%TMPR(i,j,k)  *TTREF
            cmu= BLK%VSSV(i,j,k)  *VSREF
            sc01= BLK%YSPC(i,j,k,1)
            sc02= BLK%YSPC(i,j,k,2)

            QMEN1(ii,jj,kk,1) = QMEN1(ii,jj,kk,1) + rh *rnperi
            QMEN1(ii,jj,kk,2) = QMEN1(ii,jj,kk,2) + u  *rnperi
            QMEN1(ii,jj,kk,3) = QMEN1(ii,jj,kk,3) + v  *rnperi
            QMEN1(ii,jj,kk,4) = QMEN1(ii,jj,kk,4) + w  *rnperi
            QMEN1(ii,jj,kk,5) = QMEN1(ii,jj,kk,5) + p  *rnperi
            QMEN1(ii,jj,kk,6) = QMEN1(ii,jj,kk,6) + t  *rnperi
            QMEN1(ii,jj,kk,7) = QMEN1(ii,jj,kk,7) + sc01  *rnperi
            QMEN1(ii,jj,kk,8) = QMEN1(ii,jj,kk,8) + sc02  *rnperi
          enddo
        enddo
      enddo

      if (MOD(nstep,nperiod)==0) then
        
        write(fname,900) nstep
        
        call mpisub_sbsp_out_plot3d(fname,1,QMEN1,imax,jmax,kmax,lmax)
        
        QMEN1(:,:,:,:) = 0.0d0
      end if 


    endif

    if (nstep==step_final) then
      deallocate(QMEN1)
    endif
  

  end subroutine sub_mean_Tpipe

  subroutine sub_mean(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: nstep,step_start,step_final,nperiod
    integer :: i,j,k,ii,jj,kk
    real(8) :: rh,u,v,w,p,t,cmu,sc01,sc02
    real(8) :: rnperi
    character(100) :: fname

900 format('data/mean_eta_zta_re',i7.7,'.dat')

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call steps_get_step_now(nstep)
    call steps_get_step_start(step_start)
    call steps_get_step_final(step_final)


    nperiod = 50000


    if (nstep==step_start+1) then
      allocate(QMEN2(imax,jmax,kmax,lmax))
      QMEN2(:,:,:,:) = 0.0d0
    endif


    if (step_start+1<=nstep.and.nstep<=step_final) then

      rnperi = 1.0d0/dble(nperiod)

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1

            rh = BLK%DNST(i,j,k)*RHREF
            u  = BLK%ULCT(i,j,k)*UUREF
            v  = BLK%VLCT(i,j,k)*UUREF
            w  = BLK%WLCT(i,j,k)*UUREF
            p  = BLK%PRSS(i,j,k)*PPREF
            t  = BLK%TMPR(i,j,k)*TTREF
			      sc01= BLK%YSPC(i,j,k,1)
            sc02= BLK%YSPC(i,j,k,2)

            QMEN2(ii,jj,kk,1) = QMEN2(ii,jj,kk,1) + rh *rnperi
            QMEN2(ii,jj,kk,2) = QMEN2(ii,jj,kk,2) + u  *rnperi
            QMEN2(ii,jj,kk,3) = QMEN2(ii,jj,kk,3) + v  *rnperi
            QMEN2(ii,jj,kk,4) = QMEN2(ii,jj,kk,4) + w  *rnperi
            QMEN2(ii,jj,kk,5) = QMEN2(ii,jj,kk,5) + p  *rnperi
            QMEN2(ii,jj,kk,6) = QMEN2(ii,jj,kk,6) + t  *rnperi
            QMEN2(ii,jj,kk,7) = QMEN2(ii,jj,kk,7) + sc01  *rnperi
            QMEN2(ii,jj,kk,8) = QMEN2(ii,jj,kk,8) + sc02  *rnperi
          enddo
        enddo
      enddo

    endif
    if (MOD(nstep,nperiod)==0) then
        
        write(fname,900) nstep
        
        call mpisub_sbsp_out_mean_along_eta_zta(fname,QMEN2,imax,jmax,kmax,lmax)
        
        QMEN2(:,:,:,:) = 0.0d0
    end if   

    if (nstep==step_final) then
      deallocate(QMEN2)
    endif
  
  end subroutine sub_mean

!-------------------------------------------------------------------------------
  subroutine sub_rms(BLK)
    type(block), intent(in) :: BLK

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: nstep,step_start,step_final,nperiod
    integer :: i,j,k,ii,jj,kk
    real(8) :: rh,u,v,w,p,t,cmu
    real(8) :: rnperi
    real(8) :: rh_b,u_b,v_b,w_b,p_b,t_b,cmu_b
    character(100) :: fname

900 format('data/rms_n',i7.7,'.dat')

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call steps_get_step_now(nstep)
    call steps_get_step_start(step_start)
    call steps_get_step_final(step_final)

    nperiod = step_final-(step_start+1)+1

    if (nstep==step_start+1) then

      allocate(QMEN1(imax,jmax,kmax,lmax))
      QMEN1(:,:,:,:) = 0.0d0

      allocate(QRMS(imax,jmax,kmax,lmax))
      QRMS(:,:,:,:) = 0.0d0

      call mpisub_sbsp_read_restart_from_2D_to_3D(prefix,suffix,step_rst_mean, &
           QMEN1,imax,jmax,kmax,lmax)

    endif


    if (step_start+1<=nstep.and.nstep<=step_final) then

      rnperi = 1.0d0/dble(nperiod)

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1

            rh = BLK%DNST(i,j,k)
            u  = BLK%ULCT(i,j,k)
            v  = BLK%VLCT(i,j,k)
            w  = BLK%WLCT(i,j,k)
            p  = BLK%PRSS(i,j,k)
            t  = BLK%TMPR(i,j,k)
            cmu= BLK%VSSV(i,j,k)

            rh_b = rh -QMEN1(ii,jj,kk,1)/RHREF
            u_b  = u  -QMEN1(ii,jj,kk,2)/UUREF
            v_b  = v  -QMEN1(ii,jj,kk,3)/UUREF
            w_b  = w  -QMEN1(ii,jj,kk,4)/UUREF
            p_b  = p  -QMEN1(ii,jj,kk,5)/PPREF
            t_b  = t  -QMEN1(ii,jj,kk,6)/TTREF
            cmu_b= cmu-QMEN1(ii,jj,kk,7)/VSREF

            QRMS(ii,jj,kk,1) = QRMS(ii,jj,kk,1) + rh_b* rh_b*rnperi
            QRMS(ii,jj,kk,2) = QRMS(ii,jj,kk,2) +  u_b*  u_b*rnperi
            QRMS(ii,jj,kk,3) = QRMS(ii,jj,kk,3) +  v_b*  v_b*rnperi
            QRMS(ii,jj,kk,4) = QRMS(ii,jj,kk,4) +  w_b*  w_b*rnperi
            QRMS(ii,jj,kk,5) = QRMS(ii,jj,kk,5) +  p_b*  p_b*rnperi
            QRMS(ii,jj,kk,6) = QRMS(ii,jj,kk,6) +  t_b*  t_b*rnperi
            QRMS(ii,jj,kk,7) = QRMS(ii,jj,kk,7) +cmu_b*cmu_b*rnperi
            QRMS(ii,jj,kk,8) = BLK%YNKK(i,j,k)
          enddo
        enddo
      enddo

    endif


    if (nstep==step_final) then

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1

            QRMS(ii,jj,kk,1) = dsqrt(QRMS(ii,jj,kk,1))*RHREF
            QRMS(ii,jj,kk,2) = dsqrt(QRMS(ii,jj,kk,2))*UUREF
            QRMS(ii,jj,kk,3) = dsqrt(QRMS(ii,jj,kk,3))*UUREF
            QRMS(ii,jj,kk,4) = dsqrt(QRMS(ii,jj,kk,4))*UUREF
            QRMS(ii,jj,kk,5) = dsqrt(QRMS(ii,jj,kk,5))*PPREF
            QRMS(ii,jj,kk,6) = dsqrt(QRMS(ii,jj,kk,6))*TTREF
            QRMS(ii,jj,kk,7) = dsqrt(QRMS(ii,jj,kk,7))*VSREF
            QRMS(ii,jj,kk,8) = BLK%YNKK(i,j,k)*XYREF
          enddo
        enddo
      enddo

      write(fname,900) nstep
!      call mpisub_sbsp_out_plot3d_mean_along_zta(fname,QRMS,imax,jmax,kmax,lmax)
      call mpisub_sbsp_out_mean_along_eta_zta(fname,QRMS,imax,jmax,kmax,lmax)

      deallocate(QMEN1)
      deallocate(QRMS)
    endif

  end subroutine sub_rms

  subroutine sub_reaction_rate(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax,lmax
    integer :: isplt_blk,jsplt_blk,ksplt_blk,is_blk,js_blk,ks_blk
    integer :: i,j,k,l,ii,jj,kk,im,ip,jm,jp
    real(8) :: xl,yl
    real(8) :: rh,u,v,sc02,mf
    real(8) :: mf_in1,mf_in2,mf_out
    real(8) :: u0,t0,gr,grl

    integer :: nstep,step_start

    real(8), allocatable, dimension(:,:,:,:) :: QB
    real(8), allocatable, dimension(:,:,:,:) :: QD_1
    real(8), allocatable, dimension(:,:,:,:) :: QD_2
    real(8), allocatable, dimension(:,:,:,:) :: QD_3

    900 format(9A13)
    901 format(I13,8E13.5)

    fname = 'data/reaction_rate.txt'

    lmax = 9

    u0 = UUREF*2.0d0/3.0d0
    t0 = (24.5d0*XYREF)/u0

    call mpisub_sbsp_get_myrank_world(myrank)

    call steps_get_step_now(nstep)

    call steps_get_step_start(step_start)

    if (nstep==step_start+1) then
      if (myrank==0) then
        open(69,file=fname,form='formatted',status='replace')
        write(69,900) 'STEP','MF_IN1','MF_IN2','MF_OUT','x','-ln(1-x)','t0','k','u0'
        close(69)
      endif
    endif


    if (.not.(nstep>step_start.and.mod(nstep,100)==0)) return


    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain, &
                                                   is_blk,js_blk,ks_blk)

    if (id_domain==1.and.is_blk==1) then
      allocate(QB(1,BLK%jmax,BLK%kmax,lmax))
      QB(:,:,:,:) = 0.0d0

      i = BLK%ista+1
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QB(1,jj,kk,1) = BLK%XCRD(i,j,k)
          QB(1,jj,kk,2) = BLK%YCRD(i,j,k)
          QB(1,jj,kk,3) = BLK%YSPC(i,j,k,2)
          QB(1,jj,kk,4) = BLK%DNST(i,j,k)
          QB(1,jj,kk,5) = BLK%ULCT(i,j,k)
          QB(1,jj,kk,6) = BLK%VLCT(i,j,k)
        enddo
      enddo

    endif

    if (id_domain==3.and.is_blk==isplt_blk) then
      allocate(QB(1,BLK%jmax,BLK%kmax,lmax))
      QB(:,:,:,:) = 0.0d0

      i = BLK%iend-1
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QB(1,jj,kk,1) = BLK%XCRD(i,j,k)
          QB(1,jj,kk,2) = BLK%YCRD(i,j,k)
          QB(1,jj,kk,3) = BLK%YSPC(i,j,k,2)
          QB(1,jj,kk,4) = BLK%DNST(i,j,k)
          QB(1,jj,kk,5) = BLK%ULCT(i,j,k)
          QB(1,jj,kk,6) = BLK%VLCT(i,j,k)
        enddo
      enddo

    endif

    if (id_domain==4.and.js_blk==1) then
      allocate(QB(BLK%imax,1,BLK%kmax,lmax))
      QB(:,:,:,:) = 0.0d0

      j = BLK%jsta+1
      do k=BLK%ksta,BLK%kend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QB(ii,1,kk,1) = BLK%XCRD(i,j,k)
          QB(ii,1,kk,2) = BLK%YCRD(i,j,k)
          QB(ii,1,kk,3) = BLK%YSPC(i,j,k,2)
          QB(ii,1,kk,4) = BLK%DNST(i,j,k)
          QB(ii,1,kk,5) = BLK%ULCT(i,j,k)
          QB(ii,1,kk,6) = BLK%VLCT(i,j,k)
        enddo
      enddo

    endif

    if (myrank==0) then
      call Get_DomainIJKNum_From_DomainID(1,imax,jmax,kmax)
      allocate(QD_1(1,jmax,kmax,lmax))
      call Get_DomainIJKNum_From_DomainID(3,imax,jmax,kmax)
      allocate(QD_2(1,jmax,kmax,lmax))
      call Get_DomainIJKNum_From_DomainID(4,imax,jmax,kmax)
      allocate(QD_3(imax,1,kmax,lmax))
    endif

    call Get_DomainSplit_From_DomainID(1,isplt_blk,jsplt_blk,ksplt_blk)
    call mpisub_sbsp_send_domainplane_i_to_a_rank(1,1        ,QB,0,QD_1,1,lmax,0)
    !--------------------------------------------------------------------------
    call Get_DomainSplit_From_DomainID(3,isplt_blk,jsplt_blk,ksplt_blk)
    call mpisub_sbsp_send_domainplane_i_to_a_rank(3,isplt_blk,QB,0,QD_2,1,lmax,0)
    !--------------------------------------------------------------------------
    call Get_DomainSplit_From_DomainID(4,isplt_blk,jsplt_blk,ksplt_blk)
    call mpisub_sbsp_send_domainplane_j_to_a_rank(4,1        ,QB,0,QD_3,1,lmax,0)

    if (myrank==0) then

      mf_in1 = 0.0d0
      call Get_DomainIJKNum_From_DomainID(1,imax,jmax,kmax)
      do k=1,kmax
        do j=1,jmax
          jm = max0(1   ,j-1)
          jp = min0(jmax,j+1)

          yl = ( 0.5d0*(QD_1(1,jp,k,2)+QD_1(1,j ,k,2)) &
                -0.5d0*(QD_1(1,j ,k,2)+QD_1(1,jm,k,2)) )*XYREF
          yl = dabs(yl)

          sc02=QD_1(1,j,k,3)
          rh = QD_1(1,j,k,4)*RHREF
          u  = QD_1(1,j,k,5)*UUREF
          v  = QD_1(1,j,k,6)*UUREF

          mf = rh*u*yl

          mf_in1 = mf_in1+mf
        enddo
      enddo


      mf_out = 0.0d0
      call Get_DomainIJKNum_From_DomainID(3,imax,jmax,kmax)
      do k=1,kmax
        do j=1,jmax
          jm = max0(1   ,j-1)
          jp = min0(jmax,j+1)

          yl = ( 0.5d0*(QD_2(1,jp,k,2)+QD_2(1,j ,k,2)) &
                -0.5d0*(QD_2(1,j ,k,2)+QD_2(1,jm,k,2)) )*XYREF
          yl = dabs(yl)

          sc02=QD_2(1,j,k,3)
          rh = QD_2(1,j,k,4)*RHREF
          u  = QD_2(1,j,k,5)*UUREF
          v  = QD_2(1,j,k,6)*UUREF

          mf = rh*u*yl*sc02

          mf_out = mf_out+mf
        enddo
      enddo


      mf_in2 = 0.0d0
      call Get_DomainIJKNum_From_DomainID(4,imax,jmax,kmax)
      do k=1,kmax
        do i=1,imax
          im = max0(1   ,i-1)
          ip = min0(imax,i+1)

          xl = ( 0.5d0*(QD_3(ip,1,k,1)+QD_3(i,1 ,k,1)) &
                -0.5d0*(QD_3(i,1 ,k,1)+QD_3(im,1,k,1)) )*XYREF
          xl = dabs(xl)

          sc02=QD_3(i,1,k,3)
          rh = QD_3(i,1,k,4)*RHREF
          u  = QD_3(i,1,k,5)*UUREF
          v  = QD_3(i,1,k,6)*UUREF

          mf = rh*v*xl

          mf_in2 = mf_in2+mf
        enddo
      enddo

    gr = mf_out/Min(mf_in2,0.0d0)
    grl= -log(1.0d0-gr)

    endif

    if (myrank==0) then
      open(70,file=fname,form='formatted',access='append')
      write(70,901) nstep,mf_in1,mf_in2,mf_out,gr,grl,t0,grl/t0,u0
      close(70)
    endif

    if (    id_domain==1.and.is_blk==1         &
        .or.id_domain==3.and.is_blk==isplt_blk &
        .or.id_domain==4.and.js_blk==1          ) then
      deallocate(QB)
    endif

    if (myrank==0) then
      deallocate(QD_1)
      deallocate(QD_2)
      deallocate(QD_3)
    endif
  end subroutine sub_reaction_rate

end module mdl_mean_and_rms
