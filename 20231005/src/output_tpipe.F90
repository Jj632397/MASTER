module mdl_output_tpipe
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_steps
  implicit none

contains

  subroutine output_tpipe_inflow_outflow(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax,lmax
    integer :: isplt_blk,jsplt_blk,ksplt_blk,is_blk,js_blk,ks_blk
    integer :: i,j,k,l,ii,jj,kk,im,ip,jm,jp
    real(8) :: xl,yl,zl
    real(8) :: rh,u,v,w,t,hs,ht,mf
    real(8) :: mf_in1,mf_in2,mf_out
    real(8) :: tt_in1,tt_in2,tt_out
    real(8) :: ht_in1,ht_in2,ht_out

    integer :: nstep,step_start

    real(8), allocatable, dimension(:,:,:,:) :: QB
    real(8), allocatable, dimension(:,:,:,:) :: QD_1
    real(8), allocatable, dimension(:,:,:,:) :: QD_2
    real(8), allocatable, dimension(:,:,:,:) :: QD_3

    900 format(10A13)
    901 format(I13,9E13.5)

    fname = 'data/output_inflow_outflow.txt'

    lmax = 9

    call mpisub_sbsp_get_myrank_world(myrank)

    call steps_get_step_now(nstep)

    call steps_get_step_start(step_start)

    if (nstep==step_start+1) then
      if (myrank==0) then
        open(69,file=fname,form='formatted',status='replace')
        write(69,900) 'STEP','MF_IN1','MF_IN2','MF_OUT','TS_IN1','TS_IN2','TS_OUT','HT_IN1','HT_IN2','HT_OUT'
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
          QB(1,jj,kk,3) = BLK%ZCRD(i,j,k)
          QB(1,jj,kk,4) = BLK%DNST(i,j,k)
          QB(1,jj,kk,5) = BLK%ULCT(i,j,k)
          QB(1,jj,kk,6) = BLK%VLCT(i,j,k)
          QB(1,jj,kk,7) = BLK%WLCT(i,j,k)
          QB(1,jj,kk,8) = BLK%TMPR(i,j,k)
          QB(1,jj,kk,9) = BLK%HTPS(i,j,k)
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
          QB(1,jj,kk,3) = BLK%ZCRD(i,j,k)
          QB(1,jj,kk,4) = BLK%DNST(i,j,k)
          QB(1,jj,kk,5) = BLK%ULCT(i,j,k)
          QB(1,jj,kk,6) = BLK%VLCT(i,j,k)
          QB(1,jj,kk,7) = BLK%WLCT(i,j,k)
          QB(1,jj,kk,8) = BLK%TMPR(i,j,k)
          QB(1,jj,kk,9) = BLK%HTPS(i,j,k)
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
          QB(ii,1,kk,3) = BLK%ZCRD(i,j,k)
          QB(ii,1,kk,4) = BLK%DNST(i,j,k)
          QB(ii,1,kk,5) = BLK%ULCT(i,j,k)
          QB(ii,1,kk,6) = BLK%VLCT(i,j,k)
          QB(ii,1,kk,7) = BLK%WLCT(i,j,k)
          QB(ii,1,kk,8) = BLK%TMPR(i,j,k)
          QB(ii,1,kk,9) = BLK%HTPS(i,j,k)
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
      tt_in1 = 0.0d0
      ht_in1 = 0.0d0
      call Get_DomainIJKNum_From_DomainID(1,imax,jmax,kmax)
      do k=1,kmax
        do j=1,jmax
          jm = max0(1   ,j-1)
          jp = min0(jmax,j+1)

          yl = ( 0.5d0*(QD_1(1,jp,k,2)+QD_1(1,j ,k,2)) &
                -0.5d0*(QD_1(1,j ,k,2)+QD_1(1,jm,k,2)) )*XYREF
          yl = dabs(yl)

          rh = QD_1(1,j,k,4)*RHREF
          u  = QD_1(1,j,k,5)*UUREF
          v  = QD_1(1,j,k,6)*UUREF
          w  = QD_1(1,j,k,7)*UUREF
          t  = QD_1(1,j,k,8)*TTREF
          hs = QD_1(1,j,k,9)*TTREF

          ht = hs+0.5d0*(u*u+v*v+w*w)

          mf = rh*u*yl

          mf_in1 = mf_in1+mf
          tt_in1 = tt_in1+mf*t
          ht_in1 = ht_in1+mf*ht
        enddo
      enddo

      tt_in1 = tt_in1/mf_in1


      mf_out = 0.0d0
      tt_out = 0.0d0
      ht_out = 0.0d0
      call Get_DomainIJKNum_From_DomainID(3,imax,jmax,kmax)
      do k=1,kmax
        do j=1,jmax
          jm = max0(1   ,j-1)
          jp = min0(jmax,j+1)

          yl = ( 0.5d0*(QD_2(1,jp,k,2)+QD_2(1,j ,k,2)) &
                -0.5d0*(QD_2(1,j ,k,2)+QD_2(1,jm,k,2)) )*XYREF
          yl = dabs(yl)

          rh = QD_2(1,j,k,4)*RHREF
          u  = QD_2(1,j,k,5)*UUREF
          v  = QD_2(1,j,k,6)*UUREF
          w  = QD_2(1,j,k,7)*UUREF
          t  = QD_2(1,j,k,8)*TTREF
          hs = QD_2(1,j,k,9)*TTREF

          ht = hs+0.5d0*(u*u+v*v+w*w)

          mf = rh*u*yl

          mf_out = mf_out+mf
          tt_out = tt_out+mf*t
          ht_out = ht_out+mf*ht
        enddo
      enddo

      tt_out = tt_out/mf_out


      mf_in2 = 0.0d0
      tt_in2 = 0.0d0
      ht_in2 = 0.0d0
      call Get_DomainIJKNum_From_DomainID(4,imax,jmax,kmax)
      do k=1,kmax
        do i=1,imax
          im = max0(1   ,i-1)
          ip = min0(imax,i+1)

          xl = ( 0.5d0*(QD_3(ip,1,k,1)+QD_3(i,1 ,k,1)) &
                -0.5d0*(QD_3(i,1 ,k,1)+QD_3(im,1,k,1)) )*XYREF
          xl = dabs(xl)

          rh = QD_3(i,1,k,4)*RHREF
          u  = QD_3(i,1,k,5)*UUREF
          v  = QD_3(i,1,k,6)*UUREF
          w  = QD_3(i,1,k,7)*UUREF
          t  = QD_3(i,1,k,8)*TTREF
          hs = QD_3(i,1,k,9)*TTREF

          ht = hs+0.5d0*(u*u+v*v+w*w)

          mf = rh*v*xl

          mf_in2 = mf_in2+mf
          tt_in2 = tt_in2+mf*t
          ht_in2 = ht_in2+mf*ht
        enddo
      enddo

      tt_in2 = tt_in2/mf_in2

    endif

    if (myrank==0) then
      open(70,file=fname,form='formatted',access='append')
      write(70,901) nstep,mf_in1,mf_in2,mf_out,tt_in1,tt_in2,tt_out,ht_in1,ht_in2,ht_out
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
  end subroutine output_tpipe_inflow_outflow

end module mdl_output_tpipe
