module mdl_halo
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_scmm
  use mdl_exchange_domain_halo
  implicit none

  integer, allocatable, save, dimension(:) :: IFLAG_PRIM
  integer, allocatable, save, dimension(:) :: IFLAG_GRAD
  integer, allocatable, save, dimension(:) :: IFLAG_MTRC
  integer, allocatable, save, dimension(:) :: IFLAG_RHSS
  integer, allocatable, save, dimension(:) :: IFLAG_QCS1

contains

  subroutine halo_exchange_primitive(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    integer :: iact_xsi,iact_eta,iact_zta

    integer, allocatable, dimension(:) :: IFLAG_TEMP

    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_1_M, QH_TEMP_1_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_2_M, QH_TEMP_2_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_3_M, QH_TEMP_3_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_4_M, QH_TEMP_4_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_5_M, QH_TEMP_5_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_6_M, QH_TEMP_6_P

    if (halo_width<=0) return

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainPeriodicFlag_From_DomainID(id_domain,iact_xsi,iact_eta,iact_zta)

    lmax = neqns+14

    allocate(IFLAG_TEMP(6))

    allocate(QH_TEMP_1_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_1_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_2_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_2_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_3_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_3_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_4_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_4_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_5_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_5_P(imax,jmax,halo_width,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_6_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_6_P(imax,jmax,halo_width,lmax))

    QH_TEMP_1_M(:,:,:,:) = 0.0d0
    QH_TEMP_1_P(:,:,:,:) = 0.0d0
    QH_TEMP_2_M(:,:,:,:) = 0.0d0
    QH_TEMP_2_P(:,:,:,:) = 0.0d0
    QH_TEMP_3_M(:,:,:,:) = 0.0d0
    QH_TEMP_3_P(:,:,:,:) = 0.0d0
    QH_TEMP_4_M(:,:,:,:) = 0.0d0
    QH_TEMP_4_P(:,:,:,:) = 0.0d0
    QH_TEMP_5_M(:,:,:,:) = 0.0d0
    QH_TEMP_5_P(:,:,:,:) = 0.0d0
    QH_TEMP_6_M(:,:,:,:) = 0.0d0
    QH_TEMP_6_P(:,:,:,:) = 0.0d0

    !1
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%ista+halo_width-1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_1_M(ii,jj,kk,1) = BLK%TMPR(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,2) = BLK%ULCT(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,3) = BLK%VLCT(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,4) = BLK%WLCT(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,5) = BLK%PRSS(i,j,k)
          if (neqns>=6) QH_TEMP_1_M(ii,jj,kk,6) = BLK%YSPC(i,j,k,1)
          if (neqns>=7) QH_TEMP_1_M(ii,jj,kk,7) = BLK%YSPC(i,j,k,2)
          if (neqns>=8) QH_TEMP_1_M(ii,jj,kk,8) = BLK%YSPC(i,j,k,3)
          if (neqns>=9) QH_TEMP_1_M(ii,jj,kk,9) = BLK%YSPC(i,j,k,4)
          if (neqns>=10)QH_TEMP_1_M(ii,jj,kk,10)= BLK%YSPC(i,j,k,5)
          if (neqns>=11)QH_TEMP_1_M(ii,jj,kk,11)= BLK%YSPC(i,j,k,6)
          if (neqns>=12)QH_TEMP_1_M(ii,jj,kk,12)= BLK%YSPC(i,j,k,7)
          QH_TEMP_1_M(ii,jj,kk,neqns+1) = BLK%HTPS(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+2) = BLK%DNST(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+3) = BLK%CSOS(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+4) = BLK%VSSV(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+5) = BLK%VSBV(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+6) = BLK%CKPV(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+7) = BLK%SHCP(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+8) = BLK%DTLC(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+9) = BLK%YNKK(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+10)= BLK%HFWL(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+11)= BLK%DHPR(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+12)= BLK%DRTM(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+13)= BLK%DRPR(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,neqns+14)= BLK%VPRE(i,j,k)
        enddo
      enddo
    enddo

    !2
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%iend-halo_width+1,BLK%iend
          ii = i-(BLK%iend-halo_width+1)+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_2_M(ii,jj,kk,1) = BLK%TMPR(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,2) = BLK%ULCT(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,3) = BLK%VLCT(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,4) = BLK%WLCT(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,5) = BLK%PRSS(i,j,k)
          if (neqns>=6) QH_TEMP_2_M(ii,jj,kk,6) = BLK%YSPC(i,j,k,1)
          if (neqns>=7) QH_TEMP_2_M(ii,jj,kk,7) = BLK%YSPC(i,j,k,2)
          if (neqns>=8) QH_TEMP_2_M(ii,jj,kk,8) = BLK%YSPC(i,j,k,3)
          if (neqns>=9) QH_TEMP_2_M(ii,jj,kk,9) = BLK%YSPC(i,j,k,4)
          if (neqns>=10)QH_TEMP_2_M(ii,jj,kk,10)= BLK%YSPC(i,j,k,5)
          if (neqns>=11)QH_TEMP_2_M(ii,jj,kk,11)= BLK%YSPC(i,j,k,6)
          if (neqns>=12)QH_TEMP_2_M(ii,jj,kk,12)= BLK%YSPC(i,j,k,7)
          QH_TEMP_2_M(ii,jj,kk,neqns+1) = BLK%HTPS(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+2) = BLK%DNST(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+3) = BLK%CSOS(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+4) = BLK%VSSV(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+5) = BLK%VSBV(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+6) = BLK%CKPV(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+7) = BLK%SHCP(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+8) = BLK%DTLC(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+9) = BLK%YNKK(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+10)= BLK%HFWL(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+11)= BLK%DHPR(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+12)= BLK%DRTM(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+13)= BLK%DRPR(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,neqns+14)= BLK%VPRE(i,j,k)
        enddo
      enddo
    enddo

    !3
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jsta+halo_width-1
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_3_M(ii,jj,kk,1) = BLK%TMPR(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,2) = BLK%ULCT(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,3) = BLK%VLCT(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,4) = BLK%WLCT(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,5) = BLK%PRSS(i,j,k)
          if (neqns>=6) QH_TEMP_3_M(ii,jj,kk,6) = BLK%YSPC(i,j,k,1)
          if (neqns>=7) QH_TEMP_3_M(ii,jj,kk,7) = BLK%YSPC(i,j,k,2)
          if (neqns>=8) QH_TEMP_3_M(ii,jj,kk,8) = BLK%YSPC(i,j,k,3)
          if (neqns>=9) QH_TEMP_3_M(ii,jj,kk,9) = BLK%YSPC(i,j,k,4)
          if (neqns>=10)QH_TEMP_3_M(ii,jj,kk,10)= BLK%YSPC(i,j,k,5)
          if (neqns>=11)QH_TEMP_3_M(ii,jj,kk,11)= BLK%YSPC(i,j,k,6)
          if (neqns>=12)QH_TEMP_3_M(ii,jj,kk,12)= BLK%YSPC(i,j,k,7)
          QH_TEMP_3_M(ii,jj,kk,neqns+1) = BLK%HTPS(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+2) = BLK%DNST(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+3) = BLK%CSOS(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+4) = BLK%VSSV(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+5) = BLK%VSBV(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+6) = BLK%CKPV(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+7) = BLK%SHCP(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+8) = BLK%DTLC(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+9) = BLK%YNKK(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+10)= BLK%HFWL(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+11)= BLK%DHPR(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+12)= BLK%DRTM(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+13)= BLK%DRPR(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,neqns+14)= BLK%VPRE(i,j,k)
        enddo
      enddo
    enddo

    !4
    do k=BLK%ksta,BLK%kend
      do j=BLK%jend-halo_width+1,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-(BLK%jend-halo_width+1)+1
          kk = k-BLK%ksta+1

          QH_TEMP_4_M(ii,jj,kk,1) = BLK%TMPR(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,2) = BLK%ULCT(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,3) = BLK%VLCT(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,4) = BLK%WLCT(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,5) = BLK%PRSS(i,j,k)
          if (neqns>=6) QH_TEMP_4_M(ii,jj,kk,6) = BLK%YSPC(i,j,k,1)
          if (neqns>=7) QH_TEMP_4_M(ii,jj,kk,7) = BLK%YSPC(i,j,k,2)
          if (neqns>=8) QH_TEMP_4_M(ii,jj,kk,8) = BLK%YSPC(i,j,k,3)
          if (neqns>=9) QH_TEMP_4_M(ii,jj,kk,9) = BLK%YSPC(i,j,k,4)
          if (neqns>=10)QH_TEMP_4_M(ii,jj,kk,10)= BLK%YSPC(i,j,k,5)
          if (neqns>=11)QH_TEMP_4_M(ii,jj,kk,11)= BLK%YSPC(i,j,k,6)
          if (neqns>=12)QH_TEMP_4_M(ii,jj,kk,12)= BLK%YSPC(i,j,k,7)
          QH_TEMP_4_M(ii,jj,kk,neqns+1) = BLK%HTPS(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+2) = BLK%DNST(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+3) = BLK%CSOS(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+4) = BLK%VSSV(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+5) = BLK%VSBV(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+6) = BLK%CKPV(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+7) = BLK%SHCP(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+8) = BLK%DTLC(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+9) = BLK%YNKK(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+10)= BLK%HFWL(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+11)= BLK%DHPR(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+12)= BLK%DRTM(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+13)= BLK%DRPR(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,neqns+14)= BLK%VPRE(i,j,k)
        enddo
      enddo
    enddo

    if (.not.(i2dimension>0)) then
    !5
    do k=BLK%ksta,BLK%ksta+halo_width-1
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_5_M(ii,jj,kk,1) = BLK%TMPR(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,2) = BLK%ULCT(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,3) = BLK%VLCT(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,4) = BLK%WLCT(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,5) = BLK%PRSS(i,j,k)
          if (neqns>=6) QH_TEMP_5_M(ii,jj,kk,6) = BLK%YSPC(i,j,k,1)
          if (neqns>=7) QH_TEMP_5_M(ii,jj,kk,7) = BLK%YSPC(i,j,k,2)
          if (neqns>=8) QH_TEMP_5_M(ii,jj,kk,8) = BLK%YSPC(i,j,k,3)
          if (neqns>=9) QH_TEMP_5_M(ii,jj,kk,9) = BLK%YSPC(i,j,k,4)
          if (neqns>=10)QH_TEMP_5_M(ii,jj,kk,10)= BLK%YSPC(i,j,k,5)
          if (neqns>=11)QH_TEMP_5_M(ii,jj,kk,11)= BLK%YSPC(i,j,k,6)
          if (neqns>=12)QH_TEMP_5_M(ii,jj,kk,12)= BLK%YSPC(i,j,k,7)
          QH_TEMP_5_M(ii,jj,kk,neqns+1) = BLK%HTPS(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+2) = BLK%DNST(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+3) = BLK%CSOS(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+4) = BLK%VSSV(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+5) = BLK%VSBV(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+6) = BLK%CKPV(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+7) = BLK%SHCP(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+8) = BLK%DTLC(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+9) = BLK%YNKK(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+10)= BLK%HFWL(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+11)= BLK%DHPR(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+12)= BLK%DRTM(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+13)= BLK%DRPR(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,neqns+14)= BLK%VPRE(i,j,k)
        enddo
      enddo
    enddo

    !6
    do k=BLK%kend-halo_width+1,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-(BLK%kend-halo_width+1)+1

          QH_TEMP_6_M(ii,jj,kk,1) = BLK%TMPR(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,2) = BLK%ULCT(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,3) = BLK%VLCT(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,4) = BLK%WLCT(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,5) = BLK%PRSS(i,j,k)
          if (neqns>=6) QH_TEMP_6_M(ii,jj,kk,6) = BLK%YSPC(i,j,k,1)
          if (neqns>=7) QH_TEMP_6_M(ii,jj,kk,7) = BLK%YSPC(i,j,k,2)
          if (neqns>=8) QH_TEMP_6_M(ii,jj,kk,8) = BLK%YSPC(i,j,k,3)
          if (neqns>=9) QH_TEMP_6_M(ii,jj,kk,9) = BLK%YSPC(i,j,k,4)
          if (neqns>=10)QH_TEMP_6_M(ii,jj,kk,10)= BLK%YSPC(i,j,k,5)
          if (neqns>=11)QH_TEMP_6_M(ii,jj,kk,11)= BLK%YSPC(i,j,k,6)
          if (neqns>=12)QH_TEMP_6_M(ii,jj,kk,12)= BLK%YSPC(i,j,k,7)
          QH_TEMP_6_M(ii,jj,kk,neqns+1) = BLK%HTPS(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+2) = BLK%DNST(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+3) = BLK%CSOS(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+4) = BLK%VSSV(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+5) = BLK%VSBV(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+6) = BLK%CKPV(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+7) = BLK%SHCP(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+8) = BLK%DTLC(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+9) = BLK%YNKK(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+10)= BLK%HFWL(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+11)= BLK%DHPR(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+12)= BLK%DRTM(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+13)= BLK%DRPR(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,neqns+14)= BLK%VPRE(i,j,k)
        enddo
      enddo
    enddo

    endif

    IFLAG_TEMP(:) = 0

    !exchange
    call exchange_domain_halo(QH_TEMP_1_M,QH_TEMP_1_P,QH_TEMP_2_M,QH_TEMP_2_P, &
                              QH_TEMP_3_M,QH_TEMP_3_P,QH_TEMP_4_M,QH_TEMP_4_P, &
                              QH_TEMP_5_M,QH_TEMP_5_P,QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                 QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                 QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)

    if (.not.(i2dimension>0)) then
      call mpisub_sbsp_exchange_subdomain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                   QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_xsi>0) then
      call mpisub_sbsp_exchange_domain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_eta>0) then
      call mpisub_sbsp_exchange_domain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if ((.not.(i2dimension>0)).and.(iact_zta>0)) then
      call mpisub_sbsp_exchange_domain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    IFLAG_PRIM(:) = IFLAG_TEMP(:)

    if (IFLAG_TEMP(1)==1) then
    !1
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
          ii = i-BLK%ista_wh+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          BLK%TMPR(i,j,k) = QH_TEMP_1_P(ii,jj,kk,1)
          BLK%ULCT(i,j,k) = QH_TEMP_1_P(ii,jj,kk,2)
          BLK%VLCT(i,j,k) = QH_TEMP_1_P(ii,jj,kk,3)
          BLK%WLCT(i,j,k) = QH_TEMP_1_P(ii,jj,kk,4)
          BLK%PRSS(i,j,k) = QH_TEMP_1_P(ii,jj,kk,5)
          if (neqns>=6) BLK%YSPC(i,j,k,1) = QH_TEMP_1_P(ii,jj,kk,6)
          if (neqns>=7) BLK%YSPC(i,j,k,2) = QH_TEMP_1_P(ii,jj,kk,7)
          if (neqns>=8) BLK%YSPC(i,j,k,3) = QH_TEMP_1_P(ii,jj,kk,8)
          if (neqns>=9) BLK%YSPC(i,j,k,4) = QH_TEMP_1_P(ii,jj,kk,9)
          if (neqns>=10)BLK%YSPC(i,j,k,5) = QH_TEMP_1_P(ii,jj,kk,10)
          if (neqns>=11)BLK%YSPC(i,j,k,6) = QH_TEMP_1_P(ii,jj,kk,11)
          if (neqns>=12)BLK%YSPC(i,j,k,7) = QH_TEMP_1_P(ii,jj,kk,12)
          BLK%HTPS(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+1)
          BLK%DNST(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+2)
          BLK%CSOS(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+3)
          BLK%VSSV(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+4)
          BLK%VSBV(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+5)
          BLK%CKPV(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+6)
          BLK%SHCP(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+7)
          BLK%DTLC(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+8)
          BLK%YNKK(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+9)
          BLK%HFWL(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+10)
          BLK%DHPR(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+11)
          BLK%DRTM(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+12)
          BLK%DRPR(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+13)
          BLK%VPRE(i,j,k) = QH_TEMP_1_P(ii,jj,kk,neqns+14)
        enddo
      enddo
    enddo

    endif

    if (IFLAG_TEMP(2)==1) then
    !2
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
          ii = i-(BLK%iend_wh-halo_width+1)+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          BLK%TMPR(i,j,k) = QH_TEMP_2_P(ii,jj,kk,1)
          BLK%ULCT(i,j,k) = QH_TEMP_2_P(ii,jj,kk,2)
          BLK%VLCT(i,j,k) = QH_TEMP_2_P(ii,jj,kk,3)
          BLK%WLCT(i,j,k) = QH_TEMP_2_P(ii,jj,kk,4)
          BLK%PRSS(i,j,k) = QH_TEMP_2_P(ii,jj,kk,5)
          if (neqns>=6) BLK%YSPC(i,j,k,1) = QH_TEMP_2_P(ii,jj,kk,6)
          if (neqns>=7) BLK%YSPC(i,j,k,2) = QH_TEMP_2_P(ii,jj,kk,7)
          if (neqns>=8) BLK%YSPC(i,j,k,3) = QH_TEMP_2_P(ii,jj,kk,8)
          if (neqns>=9) BLK%YSPC(i,j,k,4) = QH_TEMP_2_P(ii,jj,kk,9)
          if (neqns>=10)BLK%YSPC(i,j,k,5) = QH_TEMP_2_P(ii,jj,kk,10)
          if (neqns>=11)BLK%YSPC(i,j,k,6) = QH_TEMP_2_P(ii,jj,kk,11)
          if (neqns>=12)BLK%YSPC(i,j,k,7) = QH_TEMP_2_P(ii,jj,kk,12)
          BLK%HTPS(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+1)
          BLK%DNST(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+2)
          BLK%CSOS(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+3)
          BLK%VSSV(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+4)
          BLK%VSBV(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+5)
          BLK%CKPV(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+6)
          BLK%SHCP(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+7)
          BLK%DTLC(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+8)
          BLK%YNKK(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+9)
          BLK%HFWL(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+10)
          BLK%DHPR(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+11)
          BLK%DRTM(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+12)
          BLK%DRPR(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+13)
          BLK%VPRE(i,j,k) = QH_TEMP_2_P(ii,jj,kk,neqns+14)
        enddo
      enddo
    enddo

    endif

    if (IFLAG_TEMP(3)==1) then
    !3
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta_wh+1
          kk = k-BLK%ksta+1

          BLK%TMPR(i,j,k) = QH_TEMP_3_P(ii,jj,kk,1)
          BLK%ULCT(i,j,k) = QH_TEMP_3_P(ii,jj,kk,2)
          BLK%VLCT(i,j,k) = QH_TEMP_3_P(ii,jj,kk,3)
          BLK%WLCT(i,j,k) = QH_TEMP_3_P(ii,jj,kk,4)
          BLK%PRSS(i,j,k) = QH_TEMP_3_P(ii,jj,kk,5)
          if (neqns>=6) BLK%YSPC(i,j,k,1) = QH_TEMP_3_P(ii,jj,kk,6)
          if (neqns>=7) BLK%YSPC(i,j,k,2) = QH_TEMP_3_P(ii,jj,kk,7)
          if (neqns>=8) BLK%YSPC(i,j,k,3) = QH_TEMP_3_P(ii,jj,kk,8)
          if (neqns>=9) BLK%YSPC(i,j,k,4) = QH_TEMP_3_P(ii,jj,kk,9)
          if (neqns>=10)BLK%YSPC(i,j,k,5) = QH_TEMP_3_P(ii,jj,kk,10)
          if (neqns>=11)BLK%YSPC(i,j,k,6) = QH_TEMP_3_P(ii,jj,kk,11)
          if (neqns>=12)BLK%YSPC(i,j,k,7) = QH_TEMP_3_P(ii,jj,kk,12)
          BLK%HTPS(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+1)
          BLK%DNST(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+2)
          BLK%CSOS(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+3)
          BLK%VSSV(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+4)
          BLK%VSBV(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+5)
          BLK%CKPV(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+6)
          BLK%SHCP(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+7)
          BLK%DTLC(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+8)
          BLK%YNKK(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+9)
          BLK%HFWL(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+10)
          BLK%DHPR(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+11)
          BLK%DRTM(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+12)
          BLK%DRPR(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+13)
          BLK%VPRE(i,j,k) = QH_TEMP_3_P(ii,jj,kk,neqns+14)
        enddo
      enddo
    enddo

    endif

    if (IFLAG_TEMP(4)==1) then
    !4
    do k=BLK%ksta,BLK%kend
      do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-(BLK%jend_wh-halo_width+1)+1
          kk = k-BLK%ksta+1

          BLK%TMPR(i,j,k) = QH_TEMP_4_P(ii,jj,kk,1)
          BLK%ULCT(i,j,k) = QH_TEMP_4_P(ii,jj,kk,2)
          BLK%VLCT(i,j,k) = QH_TEMP_4_P(ii,jj,kk,3)
          BLK%WLCT(i,j,k) = QH_TEMP_4_P(ii,jj,kk,4)
          BLK%PRSS(i,j,k) = QH_TEMP_4_P(ii,jj,kk,5)
          if (neqns>=6) BLK%YSPC(i,j,k,1) = QH_TEMP_4_P(ii,jj,kk,6)
          if (neqns>=7) BLK%YSPC(i,j,k,2) = QH_TEMP_4_P(ii,jj,kk,7)
          if (neqns>=8) BLK%YSPC(i,j,k,3) = QH_TEMP_4_P(ii,jj,kk,8)
          if (neqns>=9) BLK%YSPC(i,j,k,4) = QH_TEMP_4_P(ii,jj,kk,9)
          if (neqns>=10)BLK%YSPC(i,j,k,5) = QH_TEMP_4_P(ii,jj,kk,10)
          if (neqns>=11)BLK%YSPC(i,j,k,6) = QH_TEMP_4_P(ii,jj,kk,11)
          if (neqns>=12)BLK%YSPC(i,j,k,7) = QH_TEMP_4_P(ii,jj,kk,12)
          BLK%HTPS(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+1)
          BLK%DNST(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+2)
          BLK%CSOS(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+3)
          BLK%VSSV(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+4)
          BLK%VSBV(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+5)
          BLK%CKPV(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+6)
          BLK%SHCP(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+7)
          BLK%DTLC(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+8)
          BLK%YNKK(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+9)
          BLK%HFWL(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+10)
          BLK%DHPR(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+11)
          BLK%DRTM(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+12)
          BLK%DRPR(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+13)
          BLK%VPRE(i,j,k) = QH_TEMP_4_P(ii,jj,kk,neqns+14)
        enddo
      enddo
    enddo

    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(5)==1) then
    !5
    do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta_wh+1

          BLK%TMPR(i,j,k) = QH_TEMP_5_P(ii,jj,kk,1)
          BLK%ULCT(i,j,k) = QH_TEMP_5_P(ii,jj,kk,2)
          BLK%VLCT(i,j,k) = QH_TEMP_5_P(ii,jj,kk,3)
          BLK%WLCT(i,j,k) = QH_TEMP_5_P(ii,jj,kk,4)
          BLK%PRSS(i,j,k) = QH_TEMP_5_P(ii,jj,kk,5)
          if (neqns>=6) BLK%YSPC(i,j,k,1) = QH_TEMP_5_P(ii,jj,kk,6)
          if (neqns>=7) BLK%YSPC(i,j,k,2) = QH_TEMP_5_P(ii,jj,kk,7)
          if (neqns>=8) BLK%YSPC(i,j,k,3) = QH_TEMP_5_P(ii,jj,kk,8)
          if (neqns>=9) BLK%YSPC(i,j,k,4) = QH_TEMP_5_P(ii,jj,kk,9)
          if (neqns>=10)BLK%YSPC(i,j,k,5) = QH_TEMP_5_P(ii,jj,kk,10)
          if (neqns>=11)BLK%YSPC(i,j,k,6) = QH_TEMP_5_P(ii,jj,kk,11)
          if (neqns>=12)BLK%YSPC(i,j,k,7) = QH_TEMP_5_P(ii,jj,kk,12)
          BLK%HTPS(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+1)
          BLK%DNST(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+2)
          BLK%CSOS(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+3)
          BLK%VSSV(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+4)
          BLK%VSBV(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+5)
          BLK%CKPV(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+6)
          BLK%SHCP(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+7)
          BLK%DTLC(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+8)
          BLK%YNKK(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+9)
          BLK%HFWL(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+10)
          BLK%DHPR(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+11)
          BLK%DRTM(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+12)
          BLK%DRPR(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+13)
          BLK%VPRE(i,j,k) = QH_TEMP_5_P(ii,jj,kk,neqns+14)
        enddo
      enddo
    enddo

    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(6)==1) then
    !6
    do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-(BLK%kend_wh-halo_width+1)+1

          BLK%TMPR(i,j,k) = QH_TEMP_6_P(ii,jj,kk,1)
          BLK%ULCT(i,j,k) = QH_TEMP_6_P(ii,jj,kk,2)
          BLK%VLCT(i,j,k) = QH_TEMP_6_P(ii,jj,kk,3)
          BLK%WLCT(i,j,k) = QH_TEMP_6_P(ii,jj,kk,4)
          BLK%PRSS(i,j,k) = QH_TEMP_6_P(ii,jj,kk,5)
          if (neqns>=6) BLK%YSPC(i,j,k,1) = QH_TEMP_6_P(ii,jj,kk,6)
          if (neqns>=7) BLK%YSPC(i,j,k,2) = QH_TEMP_6_P(ii,jj,kk,7)
          if (neqns>=8) BLK%YSPC(i,j,k,3) = QH_TEMP_6_P(ii,jj,kk,8)
          if (neqns>=9) BLK%YSPC(i,j,k,4) = QH_TEMP_6_P(ii,jj,kk,9)
          if (neqns>=10)BLK%YSPC(i,j,k,5) = QH_TEMP_6_P(ii,jj,kk,10)
          if (neqns>=11)BLK%YSPC(i,j,k,6) = QH_TEMP_6_P(ii,jj,kk,11)
          if (neqns>=12)BLK%YSPC(i,j,k,7) = QH_TEMP_6_P(ii,jj,kk,12)
          BLK%HTPS(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+1)
          BLK%DNST(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+2)
          BLK%CSOS(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+3)
          BLK%VSSV(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+4)
          BLK%VSBV(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+5)
          BLK%CKPV(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+6)
          BLK%SHCP(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+7)
          BLK%DTLC(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+8)
          BLK%YNKK(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+9)
          BLK%HFWL(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+10)
          BLK%DHPR(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+11)
          BLK%DRTM(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+12)
          BLK%DRPR(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+13)
          BLK%VPRE(i,j,k) = QH_TEMP_6_P(ii,jj,kk,neqns+14)
        enddo
      enddo
    enddo

    endif

    deallocate(IFLAG_TEMP)
    deallocate(QH_TEMP_1_M)
    deallocate(QH_TEMP_1_P)
    deallocate(QH_TEMP_2_M)
    deallocate(QH_TEMP_2_P)
    deallocate(QH_TEMP_3_M)
    deallocate(QH_TEMP_3_P)
    deallocate(QH_TEMP_4_M)
    deallocate(QH_TEMP_4_P)
    deallocate(QH_TEMP_5_M)
    deallocate(QH_TEMP_5_P)
    deallocate(QH_TEMP_6_M)
    deallocate(QH_TEMP_6_P)
  end subroutine halo_exchange_primitive
!-------------------------------------------------------------------------------
  subroutine halo_exchange_gradient_and_tbl_coefficient(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    integer :: iact_xsi,iact_eta,iact_zta

    integer, allocatable, dimension(:) :: IFLAG_TEMP

    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_1_M, QH_TEMP_1_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_2_M, QH_TEMP_2_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_3_M, QH_TEMP_3_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_4_M, QH_TEMP_4_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_5_M, QH_TEMP_5_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_6_M, QH_TEMP_6_P

    if (halo_width<=0) return

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainPeriodicFlag_From_DomainID(id_domain,iact_xsi,iact_eta,iact_zta)

    lmax = (neqns+1)*3+3

    allocate(IFLAG_TEMP(6))

    allocate(QH_TEMP_1_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_1_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_2_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_2_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_3_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_3_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_4_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_4_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_5_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_5_P(imax,jmax,halo_width,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_6_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_6_P(imax,jmax,halo_width,lmax))

    QH_TEMP_1_M(:,:,:,:) = 0.0d0
    QH_TEMP_1_P(:,:,:,:) = 0.0d0
    QH_TEMP_2_M(:,:,:,:) = 0.0d0
    QH_TEMP_2_P(:,:,:,:) = 0.0d0
    QH_TEMP_3_M(:,:,:,:) = 0.0d0
    QH_TEMP_3_P(:,:,:,:) = 0.0d0
    QH_TEMP_4_M(:,:,:,:) = 0.0d0
    QH_TEMP_4_P(:,:,:,:) = 0.0d0
    QH_TEMP_5_M(:,:,:,:) = 0.0d0
    QH_TEMP_5_P(:,:,:,:) = 0.0d0
    QH_TEMP_6_M(:,:,:,:) = 0.0d0
    QH_TEMP_6_P(:,:,:,:) = 0.0d0

    !1
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%ista+halo_width-1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+1)               = BLK%DTMX(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+2)               = BLK%DULX(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+3)               = BLK%DVLX(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+4)               = BLK%DWLX(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+5)               = BLK%DPRX(i,j,k)
          if (neqns>=6) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+6) = BLK%DYSX(i,j,k,1)
          if (neqns>=7) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+7) = BLK%DYSX(i,j,k,2)
          if (neqns>=8) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+8) = BLK%DYSX(i,j,k,3)
          if (neqns>=9) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+9) = BLK%DYSX(i,j,k,4)
          if (neqns>=10)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+10)= BLK%DYSX(i,j,k,5)
          if (neqns>=11)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+11)= BLK%DYSX(i,j,k,6)
          if (neqns>=12)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+12)= BLK%DYSX(i,j,k,7)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*0+neqns+1)         = BLK%DDNX(i,j,k)

          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+1)               = BLK%DTMY(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+2)               = BLK%DULY(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+3)               = BLK%DVLY(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+4)               = BLK%DWLY(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+5)               = BLK%DPRY(i,j,k)
          if (neqns>=6) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+6) = BLK%DYSY(i,j,k,1)
          if (neqns>=7) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+7) = BLK%DYSY(i,j,k,2)
          if (neqns>=8) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+8) = BLK%DYSY(i,j,k,3)
          if (neqns>=9) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+9) = BLK%DYSY(i,j,k,4)
          if (neqns>=10)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+10)= BLK%DYSY(i,j,k,5)
          if (neqns>=11)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+11)= BLK%DYSY(i,j,k,6)
          if (neqns>=12)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+12)= BLK%DYSY(i,j,k,7)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*1+neqns+1)         = BLK%DDNY(i,j,k)

          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+1)               = BLK%DTMZ(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+2)               = BLK%DULZ(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+3)               = BLK%DVLZ(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+4)               = BLK%DWLZ(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+5)               = BLK%DPRZ(i,j,k)
          if (neqns>=6) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+6) = BLK%DYSZ(i,j,k,1)
          if (neqns>=7) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+7) = BLK%DYSZ(i,j,k,2)
          if (neqns>=8) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+8) = BLK%DYSZ(i,j,k,3)
          if (neqns>=9) QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+9) = BLK%DYSZ(i,j,k,4)
          if (neqns>=10)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+10)= BLK%DYSZ(i,j,k,5)
          if (neqns>=11)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+11)= BLK%DYSZ(i,j,k,6)
          if (neqns>=12)QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+12)= BLK%DYSZ(i,j,k,7)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*2+neqns+1)         = BLK%DDNZ(i,j,k)

          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*3+1)               = BLK%VSST(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*3+2)               = BLK%VSBT(i,j,k)
          QH_TEMP_1_M(ii,jj,kk,(neqns+1)*3+3)               = BLK%CKPT(i,j,k)
        enddo
      enddo
    enddo

    !2
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%iend-halo_width+1,BLK%iend
          ii = i-(BLK%iend-halo_width+1)+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+1)               = BLK%DTMX(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+2)               = BLK%DULX(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+3)               = BLK%DVLX(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+4)               = BLK%DWLX(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+5)               = BLK%DPRX(i,j,k)
          if (neqns>=6) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+6) = BLK%DYSX(i,j,k,1)
          if (neqns>=7) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+7) = BLK%DYSX(i,j,k,2)
          if (neqns>=8) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+8) = BLK%DYSX(i,j,k,3)
          if (neqns>=9) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+9) = BLK%DYSX(i,j,k,4)
          if (neqns>=10)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+10)= BLK%DYSX(i,j,k,5)
          if (neqns>=11)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+11)= BLK%DYSX(i,j,k,6)
          if (neqns>=12)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+12)= BLK%DYSX(i,j,k,7)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*0+neqns+1)         = BLK%DDNX(i,j,k)

          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+1)               = BLK%DTMY(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+2)               = BLK%DULY(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+3)               = BLK%DVLY(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+4)               = BLK%DWLY(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+5)               = BLK%DPRY(i,j,k)
          if (neqns>=6) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+6) = BLK%DYSY(i,j,k,1)
          if (neqns>=7) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+7) = BLK%DYSY(i,j,k,2)
          if (neqns>=8) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+8) = BLK%DYSY(i,j,k,3)
          if (neqns>=9) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+9) = BLK%DYSY(i,j,k,4)
          if (neqns>=10)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+10)= BLK%DYSY(i,j,k,5)
          if (neqns>=11)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+11)= BLK%DYSY(i,j,k,6)
          if (neqns>=12)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+12)= BLK%DYSY(i,j,k,7)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*1+neqns+1)         = BLK%DDNY(i,j,k)

          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+1)               = BLK%DTMZ(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+2)               = BLK%DULZ(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+3)               = BLK%DVLZ(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+4)               = BLK%DWLZ(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+5)               = BLK%DPRZ(i,j,k)
          if (neqns>=6) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+6) = BLK%DYSZ(i,j,k,1)
          if (neqns>=7) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+7) = BLK%DYSZ(i,j,k,2)
          if (neqns>=8) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+8) = BLK%DYSZ(i,j,k,3)
          if (neqns>=9) QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+9) = BLK%DYSZ(i,j,k,4)
          if (neqns>=10)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+10)= BLK%DYSZ(i,j,k,5)
          if (neqns>=11)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+11)= BLK%DYSZ(i,j,k,6)
          if (neqns>=12)QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+12)= BLK%DYSZ(i,j,k,7)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*2+neqns+1)         = BLK%DDNZ(i,j,k)

          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*3+1)               = BLK%VSST(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*3+2)               = BLK%VSBT(i,j,k)
          QH_TEMP_2_M(ii,jj,kk,(neqns+1)*3+3)               = BLK%CKPT(i,j,k)
        enddo
      enddo
    enddo

    !3
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jsta+halo_width-1
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+1)               = BLK%DTMX(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+2)               = BLK%DULX(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+3)               = BLK%DVLX(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+4)               = BLK%DWLX(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+5)               = BLK%DPRX(i,j,k)
          if (neqns>=6) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+6) = BLK%DYSX(i,j,k,1)
          if (neqns>=7) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+7) = BLK%DYSX(i,j,k,2)
          if (neqns>=8) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+8) = BLK%DYSX(i,j,k,3)
          if (neqns>=9) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+9) = BLK%DYSX(i,j,k,4)
          if (neqns>=10)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+10)= BLK%DYSX(i,j,k,5)
          if (neqns>=11)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+11)= BLK%DYSX(i,j,k,6)
          if (neqns>=12)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+12)= BLK%DYSX(i,j,k,7)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*0+neqns+1)         = BLK%DDNX(i,j,k)

          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+1)               = BLK%DTMY(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+2)               = BLK%DULY(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+3)               = BLK%DVLY(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+4)               = BLK%DWLY(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+5)               = BLK%DPRY(i,j,k)
          if (neqns>=6) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+6) = BLK%DYSY(i,j,k,1)
          if (neqns>=7) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+7) = BLK%DYSY(i,j,k,2)
          if (neqns>=8) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+8) = BLK%DYSY(i,j,k,3)
          if (neqns>=9) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+9) = BLK%DYSY(i,j,k,4)
          if (neqns>=10)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+10)= BLK%DYSY(i,j,k,5)
          if (neqns>=11)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+11)= BLK%DYSY(i,j,k,6)
          if (neqns>=12)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+12)= BLK%DYSY(i,j,k,7)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*1+neqns+1)         = BLK%DDNY(i,j,k)

          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+1)               = BLK%DTMZ(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+2)               = BLK%DULZ(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+3)               = BLK%DVLZ(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+4)               = BLK%DWLZ(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+5)               = BLK%DPRZ(i,j,k)
          if (neqns>=6) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+6) = BLK%DYSZ(i,j,k,1)
          if (neqns>=7) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+7) = BLK%DYSZ(i,j,k,2)
          if (neqns>=8) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+8) = BLK%DYSZ(i,j,k,3)
          if (neqns>=9) QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+9) = BLK%DYSZ(i,j,k,4)
          if (neqns>=10)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+10)= BLK%DYSZ(i,j,k,5)
          if (neqns>=11)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+11)= BLK%DYSZ(i,j,k,6)
          if (neqns>=12)QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+12)= BLK%DYSZ(i,j,k,7)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*2+neqns+1)         = BLK%DDNZ(i,j,k)

          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*3+1)               = BLK%VSST(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*3+2)               = BLK%VSBT(i,j,k)
          QH_TEMP_3_M(ii,jj,kk,(neqns+1)*3+3)               = BLK%CKPT(i,j,k)
        enddo
      enddo
    enddo

    !4
    do k=BLK%ksta,BLK%kend
      do j=BLK%jend-halo_width+1,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-(BLK%jend-halo_width+1)+1
          kk = k-BLK%ksta+1

          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+1)               = BLK%DTMX(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+2)               = BLK%DULX(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+3)               = BLK%DVLX(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+4)               = BLK%DWLX(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+5)               = BLK%DPRX(i,j,k)
          if (neqns>=6) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+6) = BLK%DYSX(i,j,k,1)
          if (neqns>=7) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+7) = BLK%DYSX(i,j,k,2)
          if (neqns>=8) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+8) = BLK%DYSX(i,j,k,3)
          if (neqns>=9) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+9) = BLK%DYSX(i,j,k,4)
          if (neqns>=10)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+10)= BLK%DYSX(i,j,k,5)
          if (neqns>=11)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+11)= BLK%DYSX(i,j,k,6)
          if (neqns>=12)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+12)= BLK%DYSX(i,j,k,7)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*0+neqns+1)         = BLK%DDNX(i,j,k)

          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+1)               = BLK%DTMY(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+2)               = BLK%DULY(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+3)               = BLK%DVLY(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+4)               = BLK%DWLY(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+5)               = BLK%DPRY(i,j,k)
          if (neqns>=6) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+6) = BLK%DYSY(i,j,k,1)
          if (neqns>=7) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+7) = BLK%DYSY(i,j,k,2)
          if (neqns>=8) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+8) = BLK%DYSY(i,j,k,3)
          if (neqns>=9) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+9) = BLK%DYSY(i,j,k,4)
          if (neqns>=10)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+10)= BLK%DYSY(i,j,k,5)
          if (neqns>=11)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+11)= BLK%DYSY(i,j,k,6)
          if (neqns>=12)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+12)= BLK%DYSY(i,j,k,7)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*1+neqns+1)         = BLK%DDNY(i,j,k)

          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+1)               = BLK%DTMZ(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+2)               = BLK%DULZ(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+3)               = BLK%DVLZ(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+4)               = BLK%DWLZ(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+5)               = BLK%DPRZ(i,j,k)
          if (neqns>=6) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+6) = BLK%DYSZ(i,j,k,1)
          if (neqns>=7) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+7) = BLK%DYSZ(i,j,k,2)
          if (neqns>=8) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+8) = BLK%DYSZ(i,j,k,3)
          if (neqns>=9) QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+9) = BLK%DYSZ(i,j,k,4)
          if (neqns>=10)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+10)= BLK%DYSZ(i,j,k,5)
          if (neqns>=11)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+11)= BLK%DYSZ(i,j,k,6)
          if (neqns>=12)QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+12)= BLK%DYSZ(i,j,k,7)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*2+neqns+1)         = BLK%DDNZ(i,j,k)

          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*3+1)               = BLK%VSST(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*3+2)               = BLK%VSBT(i,j,k)
          QH_TEMP_4_M(ii,jj,kk,(neqns+1)*3+3)               = BLK%CKPT(i,j,k)
        enddo
      enddo
    enddo

    if (.not.(i2dimension>0)) then
    !5
    do k=BLK%ksta,BLK%ksta+halo_width-1
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+1)               = BLK%DTMX(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+2)               = BLK%DULX(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+3)               = BLK%DVLX(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+4)               = BLK%DWLX(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+5)               = BLK%DPRX(i,j,k)
          if (neqns>=6) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+6) = BLK%DYSX(i,j,k,1)
          if (neqns>=7) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+7) = BLK%DYSX(i,j,k,2)
          if (neqns>=8) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+8) = BLK%DYSX(i,j,k,3)
          if (neqns>=9) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+9) = BLK%DYSX(i,j,k,4)
          if (neqns>=10)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+10)= BLK%DYSX(i,j,k,5)
          if (neqns>=11)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+11)= BLK%DYSX(i,j,k,6)
          if (neqns>=12)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+12)= BLK%DYSX(i,j,k,7)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*0+neqns+1)         = BLK%DDNX(i,j,k)

          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+1)               = BLK%DTMY(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+2)               = BLK%DULY(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+3)               = BLK%DVLY(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+4)               = BLK%DWLY(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+5)               = BLK%DPRY(i,j,k)
          if (neqns>=6) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+6) = BLK%DYSY(i,j,k,1)
          if (neqns>=7) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+7) = BLK%DYSY(i,j,k,2)
          if (neqns>=8) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+8) = BLK%DYSY(i,j,k,3)
          if (neqns>=9) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+9) = BLK%DYSY(i,j,k,4)
          if (neqns>=10)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+10)= BLK%DYSY(i,j,k,5)
          if (neqns>=11)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+11)= BLK%DYSY(i,j,k,6)
          if (neqns>=12)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+12)= BLK%DYSY(i,j,k,7)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*1+neqns+1)         = BLK%DDNY(i,j,k)

          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+1)               = BLK%DTMZ(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+2)               = BLK%DULZ(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+3)               = BLK%DVLZ(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+4)               = BLK%DWLZ(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+5)               = BLK%DPRZ(i,j,k)
          if (neqns>=6) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+6) = BLK%DYSZ(i,j,k,1)
          if (neqns>=7) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+7) = BLK%DYSZ(i,j,k,2)
          if (neqns>=8) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+8) = BLK%DYSZ(i,j,k,3)
          if (neqns>=9) QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+9) = BLK%DYSZ(i,j,k,4)
          if (neqns>=10)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+10)= BLK%DYSZ(i,j,k,5)
          if (neqns>=11)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+11)= BLK%DYSZ(i,j,k,6)
          if (neqns>=12)QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+12)= BLK%DYSZ(i,j,k,7)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*2+neqns+1)         = BLK%DDNZ(i,j,k)

          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*3+1)               = BLK%VSST(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*3+2)               = BLK%VSBT(i,j,k)
          QH_TEMP_5_M(ii,jj,kk,(neqns+1)*3+3)               = BLK%CKPT(i,j,k)
        enddo
      enddo
    enddo
    endif

    if (.not.(i2dimension>0)) then
    !6
    do k=BLK%kend-halo_width+1,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-(BLK%kend-halo_width+1)+1

          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+1)               = BLK%DTMX(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+2)               = BLK%DULX(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+3)               = BLK%DVLX(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+4)               = BLK%DWLX(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+5)               = BLK%DPRX(i,j,k)
          if (neqns>=6) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+6) = BLK%DYSX(i,j,k,1)
          if (neqns>=7) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+7) = BLK%DYSX(i,j,k,2)
          if (neqns>=8) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+8) = BLK%DYSX(i,j,k,3)
          if (neqns>=9) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+9) = BLK%DYSX(i,j,k,4)
          if (neqns>=10)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+10)= BLK%DYSX(i,j,k,5)
          if (neqns>=11)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+11)= BLK%DYSX(i,j,k,6)
          if (neqns>=12)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+12)= BLK%DYSX(i,j,k,7)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*0+neqns+1)         = BLK%DDNX(i,j,k)

          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+1)               = BLK%DTMY(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+2)               = BLK%DULY(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+3)               = BLK%DVLY(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+4)               = BLK%DWLY(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+5)               = BLK%DPRY(i,j,k)
          if (neqns>=6) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+6) = BLK%DYSY(i,j,k,1)
          if (neqns>=7) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+7) = BLK%DYSY(i,j,k,2)
          if (neqns>=8) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+8) = BLK%DYSY(i,j,k,3)
          if (neqns>=9) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+9) = BLK%DYSY(i,j,k,4)
          if (neqns>=10)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+10)= BLK%DYSY(i,j,k,5)
          if (neqns>=11)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+11)= BLK%DYSY(i,j,k,6)
          if (neqns>=12)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+12)= BLK%DYSY(i,j,k,7)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*1+neqns+1)         = BLK%DDNY(i,j,k)

          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+1)               = BLK%DTMZ(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+2)               = BLK%DULZ(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+3)               = BLK%DVLZ(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+4)               = BLK%DWLZ(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+5)               = BLK%DPRZ(i,j,k)
          if (neqns>=6) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+6) = BLK%DYSZ(i,j,k,1)
          if (neqns>=7) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+7) = BLK%DYSZ(i,j,k,2)
          if (neqns>=8) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+8) = BLK%DYSZ(i,j,k,3)
          if (neqns>=9) QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+9) = BLK%DYSZ(i,j,k,4)
          if (neqns>=10)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+10)= BLK%DYSZ(i,j,k,5)
          if (neqns>=11)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+11)= BLK%DYSZ(i,j,k,6)
          if (neqns>=12)QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+12)= BLK%DYSZ(i,j,k,7)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*2+neqns+1)         = BLK%DDNZ(i,j,k)

          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*3+1)               = BLK%VSST(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*3+2)               = BLK%VSBT(i,j,k)
          QH_TEMP_6_M(ii,jj,kk,(neqns+1)*3+3)               = BLK%CKPT(i,j,k)
        enddo
      enddo
    enddo
    endif

    IFLAG_TEMP(:) = 0

    !exchange
    call exchange_domain_halo(QH_TEMP_1_M,QH_TEMP_1_P,QH_TEMP_2_M,QH_TEMP_2_P, &
                              QH_TEMP_3_M,QH_TEMP_3_P,QH_TEMP_4_M,QH_TEMP_4_P, &
                              QH_TEMP_5_M,QH_TEMP_5_P,QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP,0)


    call mpisub_sbsp_exchange_subdomain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                 QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                 QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)

    if (.not.(i2dimension>0)) then
      call mpisub_sbsp_exchange_subdomain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                   QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_xsi>0) then
      call mpisub_sbsp_exchange_domain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_eta>0) then
      call mpisub_sbsp_exchange_domain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if ((.not.(i2dimension>0)).and.(iact_zta>0)) then
      call mpisub_sbsp_exchange_domain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    IFLAG_GRAD(:) = IFLAG_TEMP(:)

    if (IFLAG_TEMP(1)==1) then
    !1
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
          ii = i-BLK%ista_wh+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          BLK%DTMX(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+1)
          BLK%DULX(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+2)
          BLK%DVLX(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+3)
          BLK%DWLX(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+4)
          BLK%DPRX(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+5)
          if (neqns>=6) BLK%DYSX(i,j,k,1) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+6)
          if (neqns>=7) BLK%DYSX(i,j,k,2) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+7)
          if (neqns>=8) BLK%DYSX(i,j,k,3) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+8)
          if (neqns>=9) BLK%DYSX(i,j,k,4) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+9)
          if (neqns>=10)BLK%DYSX(i,j,k,5) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+10)
          if (neqns>=11)BLK%DYSX(i,j,k,6) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+11)
          if (neqns>=12)BLK%DYSX(i,j,k,7) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+12)
          BLK%DDNX(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*0+neqns+1)

          BLK%DTMY(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+1)
          BLK%DULY(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+2)
          BLK%DVLY(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+3)
          BLK%DWLY(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+4)
          BLK%DPRY(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+5)
          if (neqns>=6) BLK%DYSY(i,j,k,1) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+6)
          if (neqns>=7) BLK%DYSY(i,j,k,2) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+7)
          if (neqns>=8) BLK%DYSY(i,j,k,3) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+8)
          if (neqns>=9) BLK%DYSY(i,j,k,4) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+9)
          if (neqns>=10)BLK%DYSY(i,j,k,5) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+10)
          if (neqns>=11)BLK%DYSY(i,j,k,6) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+11)
          if (neqns>=12)BLK%DYSY(i,j,k,7) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+12)
          BLK%DDNY(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*1+neqns+1)

          BLK%DTMZ(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+1)
          BLK%DULZ(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+2)
          BLK%DVLZ(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+3)
          BLK%DWLZ(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+4)
          BLK%DPRZ(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+5)
          if (neqns>=6) BLK%DYSZ(i,j,k,1) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+6)
          if (neqns>=7) BLK%DYSZ(i,j,k,2) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+7)
          if (neqns>=8) BLK%DYSZ(i,j,k,3) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+8)
          if (neqns>=9) BLK%DYSZ(i,j,k,4) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+9)
          if (neqns>=10)BLK%DYSZ(i,j,k,5) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+10)
          if (neqns>=11)BLK%DYSZ(i,j,k,6) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+11)
          if (neqns>=12)BLK%DYSZ(i,j,k,7) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+12)
          BLK%DDNZ(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*2+neqns+1)

          BLK%VSST(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*3+1)
          BLK%VSBT(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*3+2)
          BLK%CKPT(i,j,k) = QH_TEMP_1_P(ii,jj,kk,(neqns+1)*3+3)
        enddo
      enddo
    enddo

    endif

    if (IFLAG_TEMP(2)==1) then
    !1
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
          ii = i-(BLK%iend_wh-halo_width+1)+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          BLK%DTMX(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+1)
          BLK%DULX(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+2)
          BLK%DVLX(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+3)
          BLK%DWLX(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+4)
          BLK%DPRX(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+5)
          if (neqns>=6) BLK%DYSX(i,j,k,1) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+6)
          if (neqns>=7) BLK%DYSX(i,j,k,2) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+7)
          if (neqns>=8) BLK%DYSX(i,j,k,3) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+8)
          if (neqns>=9) BLK%DYSX(i,j,k,4) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+9)
          if (neqns>=10)BLK%DYSX(i,j,k,5) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+10)
          if (neqns>=11)BLK%DYSX(i,j,k,6) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+11)
          if (neqns>=12)BLK%DYSX(i,j,k,7) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+12)
          BLK%DDNX(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*0+neqns+1)

          BLK%DTMY(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+1)
          BLK%DULY(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+2)
          BLK%DVLY(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+3)
          BLK%DWLY(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+4)
          BLK%DPRY(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+5)
          if (neqns>=6) BLK%DYSY(i,j,k,1) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+6)
          if (neqns>=7) BLK%DYSY(i,j,k,2) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+7)
          if (neqns>=8) BLK%DYSY(i,j,k,3) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+8)
          if (neqns>=9) BLK%DYSY(i,j,k,4) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+9)
          if (neqns>=10)BLK%DYSY(i,j,k,5) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+10)
          if (neqns>=11)BLK%DYSY(i,j,k,6) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+11)
          if (neqns>=12)BLK%DYSY(i,j,k,7) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+12)
          BLK%DDNY(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*1+neqns+1)

          BLK%DTMZ(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+1)
          BLK%DULZ(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+2)
          BLK%DVLZ(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+3)
          BLK%DWLZ(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+4)
          BLK%DPRZ(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+5)
          if (neqns>=6) BLK%DYSZ(i,j,k,1) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+6)
          if (neqns>=7) BLK%DYSZ(i,j,k,2) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+7)
          if (neqns>=8) BLK%DYSZ(i,j,k,3) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+8)
          if (neqns>=9) BLK%DYSZ(i,j,k,4) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+9)
          if (neqns>=10)BLK%DYSZ(i,j,k,5) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+10)
          if (neqns>=11)BLK%DYSZ(i,j,k,6) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+11)
          if (neqns>=12)BLK%DYSZ(i,j,k,7) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+12)
          BLK%DDNZ(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*2+neqns+1)

          BLK%VSST(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*3+1)
          BLK%VSBT(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*3+2)
          BLK%CKPT(i,j,k) = QH_TEMP_2_P(ii,jj,kk,(neqns+1)*3+3)
        enddo
      enddo
    enddo

    endif

    if (IFLAG_TEMP(3)==1) then
    !3
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta_wh+1
          kk = k-BLK%ksta+1

          BLK%DTMX(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+1)
          BLK%DULX(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+2)
          BLK%DVLX(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+3)
          BLK%DWLX(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+4)
          BLK%DPRX(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+5)
          if (neqns>=6) BLK%DYSX(i,j,k,1) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+6)
          if (neqns>=7) BLK%DYSX(i,j,k,2) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+7)
          if (neqns>=8) BLK%DYSX(i,j,k,3) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+8)
          if (neqns>=9) BLK%DYSX(i,j,k,4) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+9)
          if (neqns>=10)BLK%DYSX(i,j,k,5) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+10)
          if (neqns>=11)BLK%DYSX(i,j,k,6) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+11)
          if (neqns>=12)BLK%DYSX(i,j,k,7) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+12)
          BLK%DDNX(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*0+neqns+1)

          BLK%DTMY(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+1)
          BLK%DULY(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+2)
          BLK%DVLY(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+3)
          BLK%DWLY(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+4)
          BLK%DPRY(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+5)
          if (neqns>=6) BLK%DYSY(i,j,k,1) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+6)
          if (neqns>=7) BLK%DYSY(i,j,k,2) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+7)
          if (neqns>=8) BLK%DYSY(i,j,k,3) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+8)
          if (neqns>=9) BLK%DYSY(i,j,k,4) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+9)
          if (neqns>=10)BLK%DYSY(i,j,k,5) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+10)
          if (neqns>=11)BLK%DYSY(i,j,k,6) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+11)
          if (neqns>=12)BLK%DYSY(i,j,k,7) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+12)
          BLK%DDNY(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*1+neqns+1)

          BLK%DTMZ(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+1)
          BLK%DULZ(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+2)
          BLK%DVLZ(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+3)
          BLK%DWLZ(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+4)
          BLK%DPRZ(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+5)
          if (neqns>=6) BLK%DYSZ(i,j,k,1) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+6)
          if (neqns>=7) BLK%DYSZ(i,j,k,2) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+7)
          if (neqns>=8) BLK%DYSZ(i,j,k,3) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+8)
          if (neqns>=9) BLK%DYSZ(i,j,k,4) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+9)
          if (neqns>=10)BLK%DYSZ(i,j,k,5) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+10)
          if (neqns>=11)BLK%DYSZ(i,j,k,6) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+11)
          if (neqns>=12)BLK%DYSZ(i,j,k,7) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+12)
          BLK%DDNZ(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*2+neqns+1)

          BLK%VSST(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*3+1)
          BLK%VSBT(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*3+2)
          BLK%CKPT(i,j,k) = QH_TEMP_3_P(ii,jj,kk,(neqns+1)*3+3)
        enddo
      enddo
    enddo

    endif

    if (IFLAG_TEMP(4)==1) then
    !4
    do k=BLK%ksta,BLK%kend
      do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-(BLK%jend_wh-halo_width+1)+1
          kk = k-BLK%ksta+1

          BLK%DTMX(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+1)
          BLK%DULX(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+2)
          BLK%DVLX(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+3)
          BLK%DWLX(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+4)
          BLK%DPRX(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+5)
          if (neqns>=6) BLK%DYSX(i,j,k,1) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+6)
          if (neqns>=7) BLK%DYSX(i,j,k,2) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+7)
          if (neqns>=8) BLK%DYSX(i,j,k,3) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+8)
          if (neqns>=9) BLK%DYSX(i,j,k,4) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+9)
          if (neqns>=10)BLK%DYSX(i,j,k,5) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+10)
          if (neqns>=11)BLK%DYSX(i,j,k,6) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+11)
          if (neqns>=12)BLK%DYSX(i,j,k,7) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+12)
          BLK%DDNX(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*0+neqns+1)

          BLK%DTMY(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+1)
          BLK%DULY(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+2)
          BLK%DVLY(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+3)
          BLK%DWLY(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+4)
          BLK%DPRY(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+5)
          if (neqns>=6) BLK%DYSY(i,j,k,1) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+6)
          if (neqns>=7) BLK%DYSY(i,j,k,2) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+7)
          if (neqns>=8) BLK%DYSY(i,j,k,3) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+8)
          if (neqns>=9) BLK%DYSY(i,j,k,4) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+9)
          if (neqns>=10)BLK%DYSY(i,j,k,5) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+10)
          if (neqns>=11)BLK%DYSY(i,j,k,6) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+11)
          if (neqns>=12)BLK%DYSY(i,j,k,7) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+12)
          BLK%DDNY(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*1+neqns+1)

          BLK%DTMZ(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+1)
          BLK%DULZ(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+2)
          BLK%DVLZ(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+3)
          BLK%DWLZ(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+4)
          BLK%DPRZ(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+5)
          if (neqns>=6) BLK%DYSZ(i,j,k,1) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+6)
          if (neqns>=7) BLK%DYSZ(i,j,k,2) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+7)
          if (neqns>=8) BLK%DYSZ(i,j,k,3) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+8)
          if (neqns>=9) BLK%DYSZ(i,j,k,4) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+9)
          if (neqns>=10)BLK%DYSZ(i,j,k,5) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+10)
          if (neqns>=11)BLK%DYSZ(i,j,k,6) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+11)
          if (neqns>=12)BLK%DYSZ(i,j,k,7) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+12)
          BLK%DDNZ(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*2+neqns+1)

          BLK%VSST(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*3+1)
          BLK%VSBT(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*3+2)
          BLK%CKPT(i,j,k) = QH_TEMP_4_P(ii,jj,kk,(neqns+1)*3+3)
        enddo
      enddo
    enddo

    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(5)==1) then
    !5
    do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta_wh+1

          BLK%DTMX(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+1)
          BLK%DULX(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+2)
          BLK%DVLX(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+3)
          BLK%DWLX(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+4)
          BLK%DPRX(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+5)
          if (neqns>=6) BLK%DYSX(i,j,k,1) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+6)
          if (neqns>=7) BLK%DYSX(i,j,k,2) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+7)
          if (neqns>=8) BLK%DYSX(i,j,k,3) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+8)
          if (neqns>=9) BLK%DYSX(i,j,k,4) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+9)
          if (neqns>=10)BLK%DYSX(i,j,k,5) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+10)
          if (neqns>=11)BLK%DYSX(i,j,k,6) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+11)
          if (neqns>=12)BLK%DYSX(i,j,k,7) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+12)
          BLK%DDNX(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*0+neqns+1)

          BLK%DTMY(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+1)
          BLK%DULY(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+2)
          BLK%DVLY(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+3)
          BLK%DWLY(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+4)
          BLK%DPRY(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+5)
          if (neqns>=6) BLK%DYSY(i,j,k,1) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+6)
          if (neqns>=7) BLK%DYSY(i,j,k,2) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+7)
          if (neqns>=8) BLK%DYSY(i,j,k,3) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+8)
          if (neqns>=9) BLK%DYSY(i,j,k,4) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+9)
          if (neqns>=10)BLK%DYSY(i,j,k,5) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+10)
          if (neqns>=11)BLK%DYSY(i,j,k,6) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+11)
          if (neqns>=12)BLK%DYSY(i,j,k,7) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+12)
          BLK%DDNY(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*1+neqns+1)

          BLK%DTMZ(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+1)
          BLK%DULZ(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+2)
          BLK%DVLZ(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+3)
          BLK%DWLZ(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+4)
          BLK%DPRZ(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+5)
          if (neqns>=6) BLK%DYSZ(i,j,k,1) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+6)
          if (neqns>=7) BLK%DYSZ(i,j,k,2) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+7)
          if (neqns>=8) BLK%DYSZ(i,j,k,3) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+8)
          if (neqns>=9) BLK%DYSZ(i,j,k,4) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+9)
          if (neqns>=10)BLK%DYSZ(i,j,k,5) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+10)
          if (neqns>=11)BLK%DYSZ(i,j,k,6) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+11)
          if (neqns>=12)BLK%DYSZ(i,j,k,7) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+12)
          BLK%DDNZ(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*2+neqns+1)

          BLK%VSST(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*3+1)
          BLK%VSBT(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*3+2)
          BLK%CKPT(i,j,k) = QH_TEMP_5_P(ii,jj,kk,(neqns+1)*3+3)
        enddo
      enddo
    enddo

    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(6)==1) then
    !6
    do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-(BLK%kend_wh-halo_width+1)+1

          BLK%DTMX(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+1)
          BLK%DULX(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+2)
          BLK%DVLX(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+3)
          BLK%DWLX(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+4)
          BLK%DPRX(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+5)
          if (neqns>=6) BLK%DYSX(i,j,k,1) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+6)
          if (neqns>=7) BLK%DYSX(i,j,k,2) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+7)
          if (neqns>=8) BLK%DYSX(i,j,k,3) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+8)
          if (neqns>=9) BLK%DYSX(i,j,k,4) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+9)
          if (neqns>=10)BLK%DYSX(i,j,k,5) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+10)
          if (neqns>=11)BLK%DYSX(i,j,k,6) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+11)
          if (neqns>=12)BLK%DYSX(i,j,k,7) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+12)
          BLK%DDNX(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*0+neqns+1)

          BLK%DTMY(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+1)
          BLK%DULY(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+2)
          BLK%DVLY(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+3)
          BLK%DWLY(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+4)
          BLK%DPRY(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+5)
          if (neqns>=6) BLK%DYSY(i,j,k,1) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+6)
          if (neqns>=7) BLK%DYSY(i,j,k,2) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+7)
          if (neqns>=8) BLK%DYSY(i,j,k,3) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+8)
          if (neqns>=9) BLK%DYSY(i,j,k,4) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+9)
          if (neqns>=10)BLK%DYSY(i,j,k,5) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+10)
          if (neqns>=11)BLK%DYSY(i,j,k,6) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+11)
          if (neqns>=12)BLK%DYSY(i,j,k,7) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+12)
          BLK%DDNY(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*1+neqns+1)

          BLK%DTMZ(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+1)
          BLK%DULZ(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+2)
          BLK%DVLZ(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+3)
          BLK%DWLZ(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+4)
          BLK%DPRZ(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+5)
          if (neqns>=6) BLK%DYSZ(i,j,k,1) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+6)
          if (neqns>=7) BLK%DYSZ(i,j,k,2) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+7)
          if (neqns>=8) BLK%DYSZ(i,j,k,3) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+8)
          if (neqns>=9) BLK%DYSZ(i,j,k,4) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+9)
          if (neqns>=10)BLK%DYSZ(i,j,k,5) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+10)
          if (neqns>=11)BLK%DYSZ(i,j,k,6) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+11)
          if (neqns>=12)BLK%DYSZ(i,j,k,7) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+12)
          BLK%DDNZ(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*2+neqns+1)

          BLK%VSST(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*3+1)
          BLK%VSBT(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*3+2)
          BLK%CKPT(i,j,k) = QH_TEMP_6_P(ii,jj,kk,(neqns+1)*3+3)
        enddo
      enddo
    enddo

    endif

    deallocate(IFLAG_TEMP)
    deallocate(QH_TEMP_1_M)
    deallocate(QH_TEMP_1_P)
    deallocate(QH_TEMP_2_M)
    deallocate(QH_TEMP_2_P)
    deallocate(QH_TEMP_3_M)
    deallocate(QH_TEMP_3_P)
    deallocate(QH_TEMP_4_M)
    deallocate(QH_TEMP_4_P)
    deallocate(QH_TEMP_5_M)
    deallocate(QH_TEMP_5_P)
    deallocate(QH_TEMP_6_M)
    deallocate(QH_TEMP_6_P)
  end subroutine halo_exchange_gradient_and_tbl_coefficient
!-------------------------------------------------------------------------------
  subroutine halo_exchange_rhs(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    integer :: iact_xsi,iact_eta,iact_zta

    integer, allocatable, dimension(:) :: IFLAG_TEMP

    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_1_M, QH_TEMP_1_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_2_M, QH_TEMP_2_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_3_M, QH_TEMP_3_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_4_M, QH_TEMP_4_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_5_M, QH_TEMP_5_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_6_M, QH_TEMP_6_P

    if (halo_width<=0) return

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainPeriodicFlag_From_DomainID(id_domain,iact_xsi,iact_eta,iact_zta)

    lmax = neqns

    allocate(IFLAG_TEMP(6))

    allocate(QH_TEMP_1_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_1_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_2_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_2_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_3_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_3_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_4_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_4_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_5_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_5_P(imax,jmax,halo_width,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_6_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_6_P(imax,jmax,halo_width,lmax))

    QH_TEMP_1_M(:,:,:,:) = 0.0d0
    QH_TEMP_1_P(:,:,:,:) = 0.0d0
    QH_TEMP_2_M(:,:,:,:) = 0.0d0
    QH_TEMP_2_P(:,:,:,:) = 0.0d0
    QH_TEMP_3_M(:,:,:,:) = 0.0d0
    QH_TEMP_3_P(:,:,:,:) = 0.0d0
    QH_TEMP_4_M(:,:,:,:) = 0.0d0
    QH_TEMP_4_P(:,:,:,:) = 0.0d0
    QH_TEMP_5_M(:,:,:,:) = 0.0d0
    QH_TEMP_5_P(:,:,:,:) = 0.0d0
    QH_TEMP_6_M(:,:,:,:) = 0.0d0
    QH_TEMP_6_P(:,:,:,:) = 0.0d0

    !1
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%ista+halo_width-1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_1_M(ii,jj,kk,l) = BLK%RHSS(i,j,k,l)/BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo

    !2
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend-halo_width+1,BLK%iend
            ii = i-(BLK%iend-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_2_M(ii,jj,kk,l) = BLK%RHSS(i,j,k,l)/BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo

    !3
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jsta+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_3_M(ii,jj,kk,l) = BLK%RHSS(i,j,k,l)/BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo

    !4
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend-halo_width+1,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend-halo_width+1)+1
            kk = k-BLK%ksta+1
            QH_TEMP_4_M(ii,jj,kk,l) = BLK%RHSS(i,j,k,l)/BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo

    if (.not.(i2dimension>0)) then
    !5
    do l=1,lmax
      do k=BLK%ksta,BLK%ksta+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_5_M(ii,jj,kk,l) = BLK%RHSS(i,j,k,l)/BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    if (.not.(i2dimension>0)) then
    !6
    do l=1,lmax
      do k=BLK%kend-halo_width+1,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend-halo_width+1)+1
            QH_TEMP_6_M(ii,jj,kk,l) = BLK%RHSS(i,j,k,l)/BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    IFLAG_TEMP(:) = 0

    !exchange
    call exchange_domain_halo(QH_TEMP_1_M,QH_TEMP_1_P,QH_TEMP_2_M,QH_TEMP_2_P, &
                              QH_TEMP_3_M,QH_TEMP_3_P,QH_TEMP_4_M,QH_TEMP_4_P, &
                              QH_TEMP_5_M,QH_TEMP_5_P,QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                 QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                 QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)

    if (.not.(i2dimension>0)) then
      call mpisub_sbsp_exchange_subdomain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                   QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_xsi>0) then
      call mpisub_sbsp_exchange_domain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_eta>0) then
      call mpisub_sbsp_exchange_domain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if ((.not.(i2dimension>0)).and.(iact_zta>0)) then
      call mpisub_sbsp_exchange_domain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    IFLAG_RHSS(:) = IFLAG_TEMP(:)

    if (IFLAG_TEMP(1)==1) then
    !1
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%RHSS(i,j,k,l) = QH_TEMP_1_P(ii,jj,kk,l)*BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(2)==1) then
    !2
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            ii = i-(BLK%iend_wh-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%RHSS(i,j,k,l) = QH_TEMP_2_P(ii,jj,kk,l)*BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(3)==1) then
    !3
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            BLK%RHSS(i,j,k,l) = QH_TEMP_3_P(ii,jj,kk,l)*BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(4)==1) then
    !4
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh-halo_width+1)+1
            kk = k-BLK%ksta+1
            BLK%RHSS(i,j,k,l) = QH_TEMP_4_P(ii,jj,kk,l)*BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(5)==1) then
    !5
    do l=1,lmax
      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            BLK%RHSS(i,j,k,l) = QH_TEMP_5_P(ii,jj,kk,l)*BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(6)==1) then
    !6
    do l=1,lmax
      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh-halo_width+1)+1
            BLK%RHSS(i,j,k,l) = QH_TEMP_6_P(ii,jj,kk,l)*BLK%RJCB(i,j,k)
          enddo
        enddo
      enddo
    enddo
    endif

    deallocate(IFLAG_TEMP)
    deallocate(QH_TEMP_1_M)
    deallocate(QH_TEMP_1_P)
    deallocate(QH_TEMP_2_M)
    deallocate(QH_TEMP_2_P)
    deallocate(QH_TEMP_3_M)
    deallocate(QH_TEMP_3_P)
    deallocate(QH_TEMP_4_M)
    deallocate(QH_TEMP_4_P)
    deallocate(QH_TEMP_5_M)
    deallocate(QH_TEMP_5_P)
    deallocate(QH_TEMP_6_M)
    deallocate(QH_TEMP_6_P)
  end subroutine halo_exchange_rhs
!-------------------------------------------------------------------------------
  subroutine halo_exchange_qcs(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    integer :: iact_xsi,iact_eta,iact_zta

    integer, allocatable, dimension(:) :: IFLAG_TEMP

    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_1_M, QH_TEMP_1_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_2_M, QH_TEMP_2_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_3_M, QH_TEMP_3_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_4_M, QH_TEMP_4_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_5_M, QH_TEMP_5_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_6_M, QH_TEMP_6_P

    if (halo_width<=0) return

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainPeriodicFlag_From_DomainID(id_domain,iact_xsi,iact_eta,iact_zta)

    lmax = neqns

    allocate(IFLAG_TEMP(6))

    allocate(QH_TEMP_1_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_1_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_2_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_2_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_3_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_3_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_4_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_4_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_5_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_5_P(imax,jmax,halo_width,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_6_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_6_P(imax,jmax,halo_width,lmax))

    QH_TEMP_1_M(:,:,:,:) = 0.0d0
    QH_TEMP_1_P(:,:,:,:) = 0.0d0
    QH_TEMP_2_M(:,:,:,:) = 0.0d0
    QH_TEMP_2_P(:,:,:,:) = 0.0d0
    QH_TEMP_3_M(:,:,:,:) = 0.0d0
    QH_TEMP_3_P(:,:,:,:) = 0.0d0
    QH_TEMP_4_M(:,:,:,:) = 0.0d0
    QH_TEMP_4_P(:,:,:,:) = 0.0d0
    QH_TEMP_5_M(:,:,:,:) = 0.0d0
    QH_TEMP_5_P(:,:,:,:) = 0.0d0
    QH_TEMP_6_M(:,:,:,:) = 0.0d0
    QH_TEMP_6_P(:,:,:,:) = 0.0d0

    !1
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%ista+halo_width-1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_1_M(ii,jj,kk,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    !2
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend-halo_width+1,BLK%iend
            ii = i-(BLK%iend-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_2_M(ii,jj,kk,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    !3
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jsta+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_3_M(ii,jj,kk,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    !4
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend-halo_width+1,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend-halo_width+1)+1
            kk = k-BLK%ksta+1
            QH_TEMP_4_M(ii,jj,kk,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    if (.not.(i2dimension>0)) then
    !5
    do l=1,lmax
      do k=BLK%ksta,BLK%ksta+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_5_M(ii,jj,kk,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (.not.(i2dimension>0)) then
    !6
    do l=1,lmax
      do k=BLK%kend-halo_width+1,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend-halo_width+1)+1
            QH_TEMP_6_M(ii,jj,kk,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
    endif

    IFLAG_TEMP(:) = 0

    !exchange
    call exchange_domain_halo(QH_TEMP_1_M,QH_TEMP_1_P,QH_TEMP_2_M,QH_TEMP_2_P, &
                              QH_TEMP_3_M,QH_TEMP_3_P,QH_TEMP_4_M,QH_TEMP_4_P, &
                              QH_TEMP_5_M,QH_TEMP_5_P,QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                 QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                 QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)

    if (.not.(i2dimension>0)) then
      call mpisub_sbsp_exchange_subdomain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                   QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_xsi>0) then
      call mpisub_sbsp_exchange_domain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_eta>0) then
      call mpisub_sbsp_exchange_domain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if ((.not.(i2dimension>0)).and.(iact_zta>0)) then
      call mpisub_sbsp_exchange_domain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    IFLAG_QCS1(:) = IFLAG_TEMP(:)

    if (IFLAG_TEMP(1)==1) then
    !1
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%QCS1(i,j,k,l) = QH_TEMP_1_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(2)==1) then
    !2
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            ii = i-(BLK%iend_wh-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%QCS1(i,j,k,l) = QH_TEMP_2_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(3)==1) then
    !3
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            BLK%QCS1(i,j,k,l) = QH_TEMP_3_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(4)==1) then
    !4
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh-halo_width+1)+1
            kk = k-BLK%ksta+1
            BLK%QCS1(i,j,k,l) = QH_TEMP_4_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(5)==1) then
    !5
    do l=1,lmax
      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            BLK%QCS1(i,j,k,l) = QH_TEMP_5_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(6)==1) then
    !6
    do l=1,lmax
      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh-halo_width+1)+1
            BLK%QCS1(i,j,k,l) = QH_TEMP_6_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    deallocate(IFLAG_TEMP)
    deallocate(QH_TEMP_1_M)
    deallocate(QH_TEMP_1_P)
    deallocate(QH_TEMP_2_M)
    deallocate(QH_TEMP_2_P)
    deallocate(QH_TEMP_3_M)
    deallocate(QH_TEMP_3_P)
    deallocate(QH_TEMP_4_M)
    deallocate(QH_TEMP_4_P)
    deallocate(QH_TEMP_5_M)
    deallocate(QH_TEMP_5_P)
    deallocate(QH_TEMP_6_M)
    deallocate(QH_TEMP_6_P)
  end subroutine halo_exchange_qcs
!-------------------------------------------------------------------------------
  subroutine halo_exchange_metric(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax
    integer :: i,j,k,ii,jj,kk

    integer :: iact_xsi,iact_eta,iact_zta

    real(8), allocatable, dimension(:,:,:,:) :: QH_XYZN_1_M, QH_XYZN_1_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_XYZN_2_M, QH_XYZN_2_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_XYZN_3_M, QH_XYZN_3_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_XYZN_4_M, QH_XYZN_4_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_XYZN_5_M, QH_XYZN_5_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_XYZN_6_M, QH_XYZN_6_P

    real(8), allocatable, dimension(:,:,:,:) :: XSIX_FACE
    real(8), allocatable, dimension(:,:,:,:) :: ETAX_FACE
    real(8), allocatable, dimension(:,:,:,:) :: ZTAX_FACE

    real(8), allocatable, dimension(:,:,:,:) :: XYZC
    real(8), allocatable, dimension(:,:,:,:) :: XYZC_FACE_ETAZTA
    real(8), allocatable, dimension(:,:,:,:) :: XYZC_FACE_ZTAXSI
    real(8), allocatable, dimension(:,:,:,:) :: XYZC_FACE_XSIETA

    real(8), allocatable, dimension(:,:,:) :: RJCB
    real(8), allocatable, dimension(:,:,:) :: GRSC

    if (halo_width<=0) return

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainPeriodicFlag_From_DomainID(id_domain,iact_xsi,iact_eta,iact_zta)

    allocate(QH_XYZN_1_M(halo_width+1,jmax+1,kmax+1,3))
    allocate(QH_XYZN_1_P(halo_width+1,jmax+1,kmax+1,3))
    allocate(QH_XYZN_2_M(halo_width+1,jmax+1,kmax+1,3))
    allocate(QH_XYZN_2_P(halo_width+1,jmax+1,kmax+1,3))
    allocate(QH_XYZN_3_M(imax+1,halo_width+1,kmax+1,3))
    allocate(QH_XYZN_3_P(imax+1,halo_width+1,kmax+1,3))
    allocate(QH_XYZN_4_M(imax+1,halo_width+1,kmax+1,3))
    allocate(QH_XYZN_4_P(imax+1,halo_width+1,kmax+1,3))
    allocate(QH_XYZN_5_M(imax+1,jmax+1,halo_width+1,3))
    allocate(QH_XYZN_5_P(imax+1,jmax+1,halo_width+1,3))
    allocate(QH_XYZN_6_M(imax+1,jmax+1,halo_width+1,3))
    allocate(QH_XYZN_6_P(imax+1,jmax+1,halo_width+1,3))

    QH_XYZN_1_M(:,:,:,:) = 0.0d0
    QH_XYZN_1_P(:,:,:,:) = 0.0d0
    QH_XYZN_2_M(:,:,:,:) = 0.0d0
    QH_XYZN_2_P(:,:,:,:) = 0.0d0
    QH_XYZN_3_M(:,:,:,:) = 0.0d0
    QH_XYZN_3_P(:,:,:,:) = 0.0d0
    QH_XYZN_4_M(:,:,:,:) = 0.0d0
    QH_XYZN_4_P(:,:,:,:) = 0.0d0
    QH_XYZN_5_M(:,:,:,:) = 0.0d0
    QH_XYZN_5_P(:,:,:,:) = 0.0d0
    QH_XYZN_6_M(:,:,:,:) = 0.0d0
    QH_XYZN_6_P(:,:,:,:) = 0.0d0

    !1
    do k=BLK%ksta,BLK%kend+1
      do j=BLK%jsta,BLK%jend+1
        do i=BLK%ista,BLK%ista+halo_width
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QH_XYZN_1_M(ii,jj,kk,1) = BLK%XCRD_NODE(i,j,k)
          QH_XYZN_1_M(ii,jj,kk,2) = BLK%YCRD_NODE(i,j,k)
          QH_XYZN_1_M(ii,jj,kk,3) = BLK%ZCRD_NODE(i,j,k)
        enddo
      enddo
    enddo

    !2
    do k=BLK%ksta,BLK%kend+1
      do j=BLK%jsta,BLK%jend+1
        do i=BLK%iend+1-halo_width,BLK%iend+1
          ii = i-(BLK%iend+1-halo_width)+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QH_XYZN_2_M(ii,jj,kk,1) = BLK%XCRD_NODE(i,j,k)
          QH_XYZN_2_M(ii,jj,kk,2) = BLK%YCRD_NODE(i,j,k)
          QH_XYZN_2_M(ii,jj,kk,3) = BLK%ZCRD_NODE(i,j,k)
        enddo
      enddo
    enddo

    !3
    do k=BLK%ksta,BLK%kend+1
      do j=BLK%jsta,BLK%jsta+halo_width
        do i=BLK%ista,BLK%iend+1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QH_XYZN_3_M(ii,jj,kk,1) = BLK%XCRD_NODE(i,j,k)
          QH_XYZN_3_M(ii,jj,kk,2) = BLK%YCRD_NODE(i,j,k)
          QH_XYZN_3_M(ii,jj,kk,3) = BLK%ZCRD_NODE(i,j,k)
        enddo
      enddo
    enddo

    !4
    do k=BLK%ksta,BLK%kend+1
      do j=BLK%jend+1-halo_width,BLK%jend+1
        do i=BLK%ista,BLK%iend+1
          ii = i-BLK%ista+1
          jj = j-(BLK%jend+1-halo_width)+1
          kk = k-BLK%ksta+1
          QH_XYZN_4_M(ii,jj,kk,1) = BLK%XCRD_NODE(i,j,k)
          QH_XYZN_4_M(ii,jj,kk,2) = BLK%YCRD_NODE(i,j,k)
          QH_XYZN_4_M(ii,jj,kk,3) = BLK%ZCRD_NODE(i,j,k)
        enddo
      enddo
    enddo

    if (.not.(i2dimension>0)) then
    !5
    do k=BLK%ksta,BLK%ksta+halo_width
      do j=BLK%jsta,BLK%jend+1
        do i=BLK%ista,BLK%iend+1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QH_XYZN_5_M(ii,jj,kk,1) = BLK%XCRD_NODE(i,j,k)
          QH_XYZN_5_M(ii,jj,kk,2) = BLK%YCRD_NODE(i,j,k)
          QH_XYZN_5_M(ii,jj,kk,3) = BLK%ZCRD_NODE(i,j,k)
        enddo
      enddo
    enddo
    endif

    if (.not.(i2dimension>0)) then
    !6
    do k=BLK%kend+1-halo_width,BLK%kend+1
      do j=BLK%jsta,BLK%jend+1
        do i=BLK%ista,BLK%iend+1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-(BLK%kend+1-halo_width)+1
          QH_XYZN_6_M(ii,jj,kk,1) = BLK%XCRD_NODE(i,j,k)
          QH_XYZN_6_M(ii,jj,kk,2) = BLK%YCRD_NODE(i,j,k)
          QH_XYZN_6_M(ii,jj,kk,3) = BLK%ZCRD_NODE(i,j,k)
        enddo
      enddo
    enddo
    endif

    IFLAG_MTRC(:) = 0

    !exchange
    call exchange_domain_halo(QH_XYZN_1_M,QH_XYZN_1_P,QH_XYZN_2_M,QH_XYZN_2_P, &
                              QH_XYZN_3_M,QH_XYZN_3_P,QH_XYZN_4_M,QH_XYZN_4_P, &
                              QH_XYZN_5_M,QH_XYZN_5_P,QH_XYZN_6_M,QH_XYZN_6_P, 3, halo_width+1, IFLAG_MTRC,1)

    call mpisub_sbsp_exchange_subdomain_halo_xsi(QH_XYZN_1_M,QH_XYZN_1_P, &
                                                 QH_XYZN_2_M,QH_XYZN_2_P, 3, halo_width+1, IFLAG_MTRC, 1)

    call mpisub_sbsp_exchange_subdomain_halo_eta(QH_XYZN_3_M,QH_XYZN_3_P, &
                                                 QH_XYZN_4_M,QH_XYZN_4_P, 3, halo_width+1, IFLAG_MTRC, 1)

    if (.not.(i2dimension>0)) then
      call mpisub_sbsp_exchange_subdomain_halo_zta(QH_XYZN_5_M,QH_XYZN_5_P, &
                                                   QH_XYZN_6_M,QH_XYZN_6_P, 3, halo_width+1, IFLAG_MTRC, 1)
    endif

    if (iact_xsi>0) then
      call mpisub_sbsp_exchange_domain_halo_xsi(QH_XYZN_1_M,QH_XYZN_1_P, &
                                                QH_XYZN_2_M,QH_XYZN_2_P, 3, halo_width+1, IFLAG_MTRC, 1)
    endif

    if (iact_eta>0) then
      call mpisub_sbsp_exchange_domain_halo_eta(QH_XYZN_3_M,QH_XYZN_3_P, &
                                                QH_XYZN_4_M,QH_XYZN_4_P, 3, halo_width+1, IFLAG_MTRC, 1)
    endif

    if ((.not.(i2dimension>0)).and.(iact_zta>0)) then
      call mpisub_sbsp_exchange_domain_halo_zta(QH_XYZN_5_M,QH_XYZN_5_P, &
                                                QH_XYZN_6_M,QH_XYZN_6_P, 3, halo_width+1, IFLAG_MTRC, 1)
    endif

    allocate(XSIX_FACE(halo_width+1,jmax,kmax,3))
    allocate(ETAX_FACE(halo_width,jmax+1,kmax,3))
    allocate(ZTAX_FACE(halo_width,jmax,kmax+1,3))
    !--------------------------------------------
    allocate(XYZC(halo_width,jmax,kmax,3))
    !---------------------------------------------------
    allocate(XYZC_FACE_ETAZTA(halo_width+1,jmax,kmax,3))
    allocate(XYZC_FACE_ZTAXSI(halo_width,jmax+1,kmax,3))
    allocate(XYZC_FACE_XSIETA(halo_width,jmax,kmax+1,3))
    !---------------------------------------------------
    allocate(RJCB(halo_width,jmax,kmax))
    allocate(GRSC(halo_width,jmax,kmax))

    !1
    if (IFLAG_MTRC(1)==1) then

      call scmm(QH_XYZN_1_P(:,:,:,1),QH_XYZN_1_P(:,:,:,2),QH_XYZN_1_P(:,:,:,3), &
                halo_width+1,jmax+1,kmax+1, &
                XSIX_FACE(:,:,:,1),XSIX_FACE(:,:,:,2),XSIX_FACE(:,:,:,3), &
                ETAX_FACE(:,:,:,1),ETAX_FACE(:,:,:,2),ETAX_FACE(:,:,:,3), &
                ZTAX_FACE(:,:,:,1),ZTAX_FACE(:,:,:,2),ZTAX_FACE(:,:,:,3), &
                RJCB)

      call node_to_cell_face(QH_XYZN_1_P(:,:,:,1),halo_width+1,jmax+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,1),XYZC_FACE_ZTAXSI(:,:,:,1),XYZC_FACE_XSIETA(:,:,:,1))

      call node_to_cell_face(QH_XYZN_1_P(:,:,:,2),halo_width+1,jmax+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,2),XYZC_FACE_ZTAXSI(:,:,:,2),XYZC_FACE_XSIETA(:,:,:,2))

      call node_to_cell_face(QH_XYZN_1_P(:,:,:,3),halo_width+1,jmax+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,3),XYZC_FACE_ZTAXSI(:,:,:,3),XYZC_FACE_XSIETA(:,:,:,3))

      call node_to_cell_center(QH_XYZN_1_P(:,:,:,1),halo_width+1,jmax+1,kmax+1,XYZC(:,:,:,1))
      call node_to_cell_center(QH_XYZN_1_P(:,:,:,2),halo_width+1,jmax+1,kmax+1,XYZC(:,:,:,2))
      call node_to_cell_center(QH_XYZN_1_P(:,:,:,3),halo_width+1,jmax+1,kmax+1,XYZC(:,:,:,3))

      call grid_scale(QH_XYZN_1_P(:,:,:,1),QH_XYZN_1_P(:,:,:,2),QH_XYZN_1_P(:,:,:,3), &
                      halo_width+1,jmax+1,kmax+1,GRSC)

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_NODE(i,j,k) = QH_XYZN_1_P(ii,jj,kk,1)
            BLK%YCRD_NODE(i,j,k) = QH_XYZN_1_P(ii,jj,kk,2)
            BLK%ZCRD_NODE(i,j,k) = QH_XYZN_1_P(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,1)
            BLK%YCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,2)
            BLK%ZCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,3)
            BLK%DXSX_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,1)
            BLK%DXSY_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,2)
            BLK%DXSZ_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DXSZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,1)
            BLK%YCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,2)
            BLK%ZCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,3)
            BLK%DETX_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,1)
            BLK%DETY_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,2)
            BLK%DETZ_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DETZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,1)
            BLK%YCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,2)
            BLK%ZCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,3)
            BLK%DZTX_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,1)
            BLK%DZTY_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,2)
            BLK%DZTZ_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DZTX_FACE(i,j,k) = 0.0d0
              BLK%DZTY_FACE(i,j,k) = 0.0d0
              BLK%DZTZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%DXSX(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,1)+XSIX_FACE(ii,jj,kk,1))
            BLK%DXSY(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,2)+XSIX_FACE(ii,jj,kk,2))
            BLK%DXSZ(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,3)+XSIX_FACE(ii,jj,kk,3))
            BLK%DETX(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,1)+ETAX_FACE(ii,jj,kk,1))
            BLK%DETY(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,2)+ETAX_FACE(ii,jj,kk,2))
            BLK%DETZ(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,3)+ETAX_FACE(ii,jj,kk,3))
            BLK%DZTX(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,1)+ZTAX_FACE(ii,jj,kk,1))
            BLK%DZTY(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,2)+ZTAX_FACE(ii,jj,kk,2))
            BLK%DZTZ(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,3)+ZTAX_FACE(ii,jj,kk,3))
            BLK%XCRD(i,j,k) = XYZC(ii,jj,kk,1)
            BLK%YCRD(i,j,k) = XYZC(ii,jj,kk,2)
            BLK%ZCRD(i,j,k) = XYZC(ii,jj,kk,3)
            BLK%RJCB(i,j,k) = RJCB(ii,jj,kk)
            BLK%GRSC(i,j,k) = GRSC(ii,jj,kk)
          enddo
        enddo
      enddo

    endif

    !2
    if (IFLAG_MTRC(2)==1) then

      call scmm(QH_XYZN_2_P(:,:,:,1),QH_XYZN_2_P(:,:,:,2),QH_XYZN_2_P(:,:,:,3), &
                halo_width+1,jmax+1,kmax+1, &
                XSIX_FACE(:,:,:,1),XSIX_FACE(:,:,:,2),XSIX_FACE(:,:,:,3), &
                ETAX_FACE(:,:,:,1),ETAX_FACE(:,:,:,2),ETAX_FACE(:,:,:,3), &
                ZTAX_FACE(:,:,:,1),ZTAX_FACE(:,:,:,2),ZTAX_FACE(:,:,:,3), &
                RJCB)


      call node_to_cell_face(QH_XYZN_2_P(:,:,:,1),halo_width+1,jmax+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,1),XYZC_FACE_ZTAXSI(:,:,:,1),XYZC_FACE_XSIETA(:,:,:,1))

      call node_to_cell_face(QH_XYZN_2_P(:,:,:,2),halo_width+1,jmax+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,2),XYZC_FACE_ZTAXSI(:,:,:,2),XYZC_FACE_XSIETA(:,:,:,2))

      call node_to_cell_face(QH_XYZN_2_P(:,:,:,3),halo_width+1,jmax+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,3),XYZC_FACE_ZTAXSI(:,:,:,3),XYZC_FACE_XSIETA(:,:,:,3))

      call node_to_cell_center(QH_XYZN_2_P(:,:,:,1),halo_width+1,jmax+1,kmax+1,XYZC(:,:,:,1))
      call node_to_cell_center(QH_XYZN_2_P(:,:,:,2),halo_width+1,jmax+1,kmax+1,XYZC(:,:,:,2))
      call node_to_cell_center(QH_XYZN_2_P(:,:,:,3),halo_width+1,jmax+1,kmax+1,XYZC(:,:,:,3))

      call grid_scale(QH_XYZN_2_P(:,:,:,1),QH_XYZN_2_P(:,:,:,2),QH_XYZN_2_P(:,:,:,3), &
                      halo_width+1,jmax+1,kmax+1,GRSC)

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%iend_wh+1-halo_width+1,BLK%iend_wh+1
            ii = i-(BLK%iend_wh+1-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_NODE(i,j,k) = QH_XYZN_2_P(ii+1,jj,kk,1)
            BLK%YCRD_NODE(i,j,k) = QH_XYZN_2_P(ii+1,jj,kk,2)
            BLK%ZCRD_NODE(i,j,k) = QH_XYZN_2_P(ii+1,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh+1-halo_width+1,BLK%iend_wh+1
            ii = i-(BLK%iend_wh+1-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii+1,jj,kk,1)
            BLK%YCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii+1,jj,kk,2)
            BLK%ZCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii+1,jj,kk,3)
            BLK%DXSX_FACE(i,j,k) = XSIX_FACE(ii+1,jj,kk,1)
            BLK%DXSY_FACE(i,j,k) = XSIX_FACE(ii+1,jj,kk,2)
            BLK%DXSZ_FACE(i,j,k) = XSIX_FACE(ii+1,jj,kk,3)
            if (i2dimension>0) then
              BLK%DXSZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            ii = i-(BLK%iend_wh-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,1)
            BLK%YCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,2)
            BLK%ZCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,3)
            BLK%DETX_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,1)
            BLK%DETY_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,2)
            BLK%DETZ_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DETZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            ii = i-(BLK%iend_wh-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,1)
            BLK%YCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,2)
            BLK%ZCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,3)
            BLK%DZTX_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,1)
            BLK%DZTY_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,2)
            BLK%DZTZ_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DZTX_FACE(i,j,k) = 0.0d0
              BLK%DZTY_FACE(i,j,k) = 0.0d0
              BLK%DZTZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            ii = i-(BLK%iend_wh-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%DXSX(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,1)+XSIX_FACE(ii,jj,kk,1))
            BLK%DXSY(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,2)+XSIX_FACE(ii,jj,kk,2))
            BLK%DXSZ(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,3)+XSIX_FACE(ii,jj,kk,3))
            BLK%DETX(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,1)+ETAX_FACE(ii,jj,kk,1))
            BLK%DETY(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,2)+ETAX_FACE(ii,jj,kk,2))
            BLK%DETZ(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,3)+ETAX_FACE(ii,jj,kk,3))
            BLK%DZTX(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,1)+ZTAX_FACE(ii,jj,kk,1))
            BLK%DZTY(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,2)+ZTAX_FACE(ii,jj,kk,2))
            BLK%DZTZ(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,3)+ZTAX_FACE(ii,jj,kk,3))
            BLK%XCRD(i,j,k) = XYZC(ii,jj,kk,1)
            BLK%YCRD(i,j,k) = XYZC(ii,jj,kk,2)
            BLK%ZCRD(i,j,k) = XYZC(ii,jj,kk,3)
            BLK%RJCB(i,j,k) = RJCB(ii,jj,kk)
            BLK%GRSC(i,j,k) = GRSC(ii,jj,kk)
          enddo
        enddo
      enddo

    endif

    deallocate(XSIX_FACE)
    deallocate(ETAX_FACE)
    deallocate(ZTAX_FACE)
    deallocate(XYZC)
    deallocate(XYZC_FACE_ETAZTA)
    deallocate(XYZC_FACE_ZTAXSI)
    deallocate(XYZC_FACE_XSIETA)
    deallocate(RJCB)
    deallocate(GRSC)


    allocate(XSIX_FACE(imax+1,halo_width,kmax,3))
    allocate(ETAX_FACE(imax,halo_width+1,kmax,3))
    allocate(ZTAX_FACE(imax,halo_width,kmax+1,3))
    !--------------------------------------------
    allocate(XYZC(imax,halo_width,kmax,3))
    !---------------------------------------------------
    allocate(XYZC_FACE_ETAZTA(imax+1,halo_width,kmax,3))
    allocate(XYZC_FACE_ZTAXSI(imax,halo_width+1,kmax,3))
    allocate(XYZC_FACE_XSIETA(imax,halo_width,kmax+1,3))
    !---------------------------------------------------
    allocate(RJCB(imax,halo_width,kmax))
    allocate(GRSC(imax,halo_width,kmax))

    !3
    if (IFLAG_MTRC(3)==1) then

      call scmm(QH_XYZN_3_P(:,:,:,1),QH_XYZN_3_P(:,:,:,2),QH_XYZN_3_P(:,:,:,3), &
                imax+1,halo_width+1,kmax+1, &
                XSIX_FACE(:,:,:,1),XSIX_FACE(:,:,:,2),XSIX_FACE(:,:,:,3), &
                ETAX_FACE(:,:,:,1),ETAX_FACE(:,:,:,2),ETAX_FACE(:,:,:,3), &
                ZTAX_FACE(:,:,:,1),ZTAX_FACE(:,:,:,2),ZTAX_FACE(:,:,:,3), &
                RJCB)

      call node_to_cell_face(QH_XYZN_3_P(:,:,:,1),imax+1,halo_width+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,1),XYZC_FACE_ZTAXSI(:,:,:,1),XYZC_FACE_XSIETA(:,:,:,1))

      call node_to_cell_face(QH_XYZN_3_P(:,:,:,2),imax+1,halo_width+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,2),XYZC_FACE_ZTAXSI(:,:,:,2),XYZC_FACE_XSIETA(:,:,:,2))

      call node_to_cell_face(QH_XYZN_3_P(:,:,:,3),imax+1,halo_width+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,3),XYZC_FACE_ZTAXSI(:,:,:,3),XYZC_FACE_XSIETA(:,:,:,3))

      call node_to_cell_center(QH_XYZN_3_P(:,:,:,1),imax+1,halo_width+1,kmax+1,XYZC(:,:,:,1))
      call node_to_cell_center(QH_XYZN_3_P(:,:,:,2),imax+1,halo_width+1,kmax+1,XYZC(:,:,:,2))
      call node_to_cell_center(QH_XYZN_3_P(:,:,:,3),imax+1,halo_width+1,kmax+1,XYZC(:,:,:,3))

      call grid_scale(QH_XYZN_3_P(:,:,:,1),QH_XYZN_3_P(:,:,:,2),QH_XYZN_3_P(:,:,:,3), &
                      imax+1,halo_width+1,kmax+1,GRSC)

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            BLK%XCRD_NODE(i,j,k) = QH_XYZN_3_P(ii,jj,kk,1)
            BLK%YCRD_NODE(i,j,k) = QH_XYZN_3_P(ii,jj,kk,2)
            BLK%ZCRD_NODE(i,j,k) = QH_XYZN_3_P(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,1)
            BLK%YCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,2)
            BLK%ZCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,3)
            BLK%DXSX_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,1)
            BLK%DXSY_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,2)
            BLK%DXSZ_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DXSZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,1)
            BLK%YCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,2)
            BLK%ZCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,3)
            BLK%DETX_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,1)
            BLK%DETY_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,2)
            BLK%DETZ_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DETZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,1)
            BLK%YCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,2)
            BLK%ZCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,3)
            BLK%DZTX_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,1)
            BLK%DZTY_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,2)
            BLK%DZTZ_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DZTX_FACE(i,j,k) = 0.0d0
              BLK%DZTY_FACE(i,j,k) = 0.0d0
              BLK%DZTZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            BLK%DXSX(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,1)+XSIX_FACE(ii,jj,kk,1))
            BLK%DXSY(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,2)+XSIX_FACE(ii,jj,kk,2))
            BLK%DXSZ(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,3)+XSIX_FACE(ii,jj,kk,3))
            BLK%DETX(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,1)+ETAX_FACE(ii,jj,kk,1))
            BLK%DETY(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,2)+ETAX_FACE(ii,jj,kk,2))
            BLK%DETZ(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,3)+ETAX_FACE(ii,jj,kk,3))
            BLK%DZTX(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,1)+ZTAX_FACE(ii,jj,kk,1))
            BLK%DZTY(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,2)+ZTAX_FACE(ii,jj,kk,2))
            BLK%DZTZ(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,3)+ZTAX_FACE(ii,jj,kk,3))
            BLK%XCRD(i,j,k) = XYZC(ii,jj,kk,1)
            BLK%YCRD(i,j,k) = XYZC(ii,jj,kk,2)
            BLK%ZCRD(i,j,k) = XYZC(ii,jj,kk,3)
            BLK%RJCB(i,j,k) = RJCB(ii,jj,kk)
            BLK%GRSC(i,j,k) = GRSC(ii,jj,kk)
          enddo
        enddo
      enddo

    endif

    !4
    if (IFLAG_MTRC(4)==1) then

      call scmm(QH_XYZN_4_P(:,:,:,1),QH_XYZN_4_P(:,:,:,2),QH_XYZN_4_P(:,:,:,3), &
                imax+1,halo_width+1,kmax+1, &
                XSIX_FACE(:,:,:,1),XSIX_FACE(:,:,:,2),XSIX_FACE(:,:,:,3), &
                ETAX_FACE(:,:,:,1),ETAX_FACE(:,:,:,2),ETAX_FACE(:,:,:,3), &
                ZTAX_FACE(:,:,:,1),ZTAX_FACE(:,:,:,2),ZTAX_FACE(:,:,:,3), &
                RJCB)


      call node_to_cell_face(QH_XYZN_4_P(:,:,:,1),imax+1,halo_width+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,1),XYZC_FACE_ZTAXSI(:,:,:,1),XYZC_FACE_XSIETA(:,:,:,1))

      call node_to_cell_face(QH_XYZN_4_P(:,:,:,2),imax+1,halo_width+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,2),XYZC_FACE_ZTAXSI(:,:,:,2),XYZC_FACE_XSIETA(:,:,:,2))

      call node_to_cell_face(QH_XYZN_4_P(:,:,:,3),imax+1,halo_width+1,kmax+1, &
           XYZC_FACE_ETAZTA(:,:,:,3),XYZC_FACE_ZTAXSI(:,:,:,3),XYZC_FACE_XSIETA(:,:,:,3))

      call node_to_cell_center(QH_XYZN_4_P(:,:,:,1),imax+1,halo_width+1,kmax+1,XYZC(:,:,:,1))
      call node_to_cell_center(QH_XYZN_4_P(:,:,:,2),imax+1,halo_width+1,kmax+1,XYZC(:,:,:,2))
      call node_to_cell_center(QH_XYZN_4_P(:,:,:,3),imax+1,halo_width+1,kmax+1,XYZC(:,:,:,3))

      call grid_scale(QH_XYZN_4_P(:,:,:,1),QH_XYZN_4_P(:,:,:,2),QH_XYZN_4_P(:,:,:,3), &
                      imax+1,halo_width+1,kmax+1,GRSC)

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jend_wh+1-halo_width+1,BLK%jend_wh+1
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh+1-halo_width+1)+1
            kk = k-BLK%ksta+1
            BLK%XCRD_NODE(i,j,k) = QH_XYZN_4_P(ii,jj+1,kk,1)
            BLK%YCRD_NODE(i,j,k) = QH_XYZN_4_P(ii,jj+1,kk,2)
            BLK%ZCRD_NODE(i,j,k) = QH_XYZN_4_P(ii,jj+1,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh-halo_width+1)+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,1)
            BLK%YCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,2)
            BLK%ZCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,3)
            BLK%DXSX_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,1)
            BLK%DXSY_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,2)
            BLK%DXSZ_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DXSZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh+1-halo_width+1,BLK%jend_wh+1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh+1-halo_width+1)+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj+1,kk,1)
            BLK%YCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj+1,kk,2)
            BLK%ZCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj+1,kk,3)
            BLK%DETX_FACE(i,j,k) = ETAX_FACE(ii,jj+1,kk,1)
            BLK%DETY_FACE(i,j,k) = ETAX_FACE(ii,jj+1,kk,2)
            BLK%DETZ_FACE(i,j,k) = ETAX_FACE(ii,jj+1,kk,3)
            if (i2dimension>0) then
              BLK%DETZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend+1
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh-halo_width+1)+1
            kk = k-BLK%ksta+1
            BLK%XCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,1)
            BLK%YCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,2)
            BLK%ZCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,3)
            BLK%DZTX_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,1)
            BLK%DZTY_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,2)
            BLK%DZTZ_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,3)
            if (i2dimension>0) then
              BLK%DZTX_FACE(i,j,k) = 0.0d0
              BLK%DZTY_FACE(i,j,k) = 0.0d0
              BLK%DZTZ_FACE(i,j,k) = 0.0d0
            endif
          enddo
        enddo
      enddo

      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh-halo_width+1)+1
            kk = k-BLK%ksta+1
            BLK%DXSX(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,1)+XSIX_FACE(ii,jj,kk,1))
            BLK%DXSY(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,2)+XSIX_FACE(ii,jj,kk,2))
            BLK%DXSZ(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,3)+XSIX_FACE(ii,jj,kk,3))
            BLK%DETX(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,1)+ETAX_FACE(ii,jj,kk,1))
            BLK%DETY(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,2)+ETAX_FACE(ii,jj,kk,2))
            BLK%DETZ(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,3)+ETAX_FACE(ii,jj,kk,3))
            BLK%DZTX(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,1)+ZTAX_FACE(ii,jj,kk,1))
            BLK%DZTY(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,2)+ZTAX_FACE(ii,jj,kk,2))
            BLK%DZTZ(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,3)+ZTAX_FACE(ii,jj,kk,3))
            BLK%XCRD(i,j,k) = XYZC(ii,jj,kk,1)
            BLK%YCRD(i,j,k) = XYZC(ii,jj,kk,2)
            BLK%ZCRD(i,j,k) = XYZC(ii,jj,kk,3)
            BLK%RJCB(i,j,k) = RJCB(ii,jj,kk)
            BLK%GRSC(i,j,k) = GRSC(ii,jj,kk)
          enddo
        enddo
      enddo

    endif

    deallocate(XSIX_FACE)
    deallocate(ETAX_FACE)
    deallocate(ZTAX_FACE)
    deallocate(XYZC)
    deallocate(XYZC_FACE_ETAZTA)
    deallocate(XYZC_FACE_ZTAXSI)
    deallocate(XYZC_FACE_XSIETA)
    deallocate(RJCB)
    deallocate(GRSC)


    allocate(XSIX_FACE(imax+1,jmax,halo_width,3))
    allocate(ETAX_FACE(imax,jmax+1,halo_width,3))
    allocate(ZTAX_FACE(imax,jmax,halo_width+1,3))
    !--------------------------------------------
    allocate(XYZC(imax,jmax,halo_width,3))
    !---------------------------------------------------
    allocate(XYZC_FACE_ETAZTA(imax+1,jmax,halo_width,3))
    allocate(XYZC_FACE_ZTAXSI(imax,jmax+1,halo_width,3))
    allocate(XYZC_FACE_XSIETA(imax,jmax,halo_width+1,3))
    !---------------------------------------------------
    allocate(RJCB(imax,jmax,halo_width))
    allocate(GRSC(imax,jmax,halo_width))

    !5
    if ((.not.(i2dimension>0)).and.(IFLAG_MTRC(5)==1)) then

      call scmm(QH_XYZN_5_P(:,:,:,1),QH_XYZN_5_P(:,:,:,2),QH_XYZN_5_P(:,:,:,3), &
                imax+1,jmax+1,halo_width+1, &
                XSIX_FACE(:,:,:,1),XSIX_FACE(:,:,:,2),XSIX_FACE(:,:,:,3), &
                ETAX_FACE(:,:,:,1),ETAX_FACE(:,:,:,2),ETAX_FACE(:,:,:,3), &
                ZTAX_FACE(:,:,:,1),ZTAX_FACE(:,:,:,2),ZTAX_FACE(:,:,:,3), &
                RJCB)

      call node_to_cell_face(QH_XYZN_5_P(:,:,:,1),imax+1,jmax+1,halo_width+1, &
           XYZC_FACE_ETAZTA(:,:,:,1),XYZC_FACE_ZTAXSI(:,:,:,1),XYZC_FACE_XSIETA(:,:,:,1))

      call node_to_cell_face(QH_XYZN_5_P(:,:,:,2),imax+1,jmax+1,halo_width+1, &
           XYZC_FACE_ETAZTA(:,:,:,2),XYZC_FACE_ZTAXSI(:,:,:,2),XYZC_FACE_XSIETA(:,:,:,2))

      call node_to_cell_face(QH_XYZN_5_P(:,:,:,3),imax+1,jmax+1,halo_width+1, &
           XYZC_FACE_ETAZTA(:,:,:,3),XYZC_FACE_ZTAXSI(:,:,:,3),XYZC_FACE_XSIETA(:,:,:,3))

      call node_to_cell_center(QH_XYZN_5_P(:,:,:,1),imax+1,jmax+1,halo_width+1,XYZC(:,:,:,1))
      call node_to_cell_center(QH_XYZN_5_P(:,:,:,2),imax+1,jmax+1,halo_width+1,XYZC(:,:,:,2))
      call node_to_cell_center(QH_XYZN_5_P(:,:,:,3),imax+1,jmax+1,halo_width+1,XYZC(:,:,:,3))

      call grid_scale(QH_XYZN_5_P(:,:,:,1),QH_XYZN_5_P(:,:,:,2),QH_XYZN_5_P(:,:,:,3), &
                      imax+1,jmax+1,halo_width+1,GRSC)

      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            BLK%XCRD_NODE(i,j,k) = QH_XYZN_5_P(ii,jj,kk,1)
            BLK%YCRD_NODE(i,j,k) = QH_XYZN_5_P(ii,jj,kk,2)
            BLK%ZCRD_NODE(i,j,k) = QH_XYZN_5_P(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            BLK%XCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,1)
            BLK%YCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,2)
            BLK%ZCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,3)
            BLK%DXSX_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,1)
            BLK%DXSY_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,2)
            BLK%DXSZ_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            BLK%XCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,1)
            BLK%YCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,2)
            BLK%ZCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,3)
            BLK%DETX_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,1)
            BLK%DETY_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,2)
            BLK%DETZ_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            BLK%XCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,1)
            BLK%YCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,2)
            BLK%ZCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk,3)
            BLK%DZTX_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,1)
            BLK%DZTY_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,2)
            BLK%DZTZ_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            BLK%DXSX(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,1)+XSIX_FACE(ii,jj,kk,1))
            BLK%DXSY(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,2)+XSIX_FACE(ii,jj,kk,2))
            BLK%DXSZ(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,3)+XSIX_FACE(ii,jj,kk,3))
            BLK%DETX(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,1)+ETAX_FACE(ii,jj,kk,1))
            BLK%DETY(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,2)+ETAX_FACE(ii,jj,kk,2))
            BLK%DETZ(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,3)+ETAX_FACE(ii,jj,kk,3))
            BLK%DZTX(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,1)+ZTAX_FACE(ii,jj,kk,1))
            BLK%DZTY(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,2)+ZTAX_FACE(ii,jj,kk,2))
            BLK%DZTZ(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,3)+ZTAX_FACE(ii,jj,kk,3))
            BLK%XCRD(i,j,k) = XYZC(ii,jj,kk,1)
            BLK%YCRD(i,j,k) = XYZC(ii,jj,kk,2)
            BLK%ZCRD(i,j,k) = XYZC(ii,jj,kk,3)
            BLK%RJCB(i,j,k) = RJCB(ii,jj,kk)
            BLK%GRSC(i,j,k) = GRSC(ii,jj,kk)
          enddo
        enddo
      enddo

    endif

    !6
    if ((.not.(i2dimension>0)).and.(IFLAG_MTRC(6)==1)) then

      call scmm(QH_XYZN_6_P(:,:,:,1),QH_XYZN_6_P(:,:,:,2),QH_XYZN_6_P(:,:,:,3), &
                imax+1,jmax+1,halo_width+1, &
                XSIX_FACE(:,:,:,1),XSIX_FACE(:,:,:,2),XSIX_FACE(:,:,:,3), &
                ETAX_FACE(:,:,:,1),ETAX_FACE(:,:,:,2),ETAX_FACE(:,:,:,3), &
                ZTAX_FACE(:,:,:,1),ZTAX_FACE(:,:,:,2),ZTAX_FACE(:,:,:,3), &
                RJCB)


      call node_to_cell_face(QH_XYZN_6_P(:,:,:,1),imax+1,jmax+1,halo_width+1, &
           XYZC_FACE_ETAZTA(:,:,:,1),XYZC_FACE_ZTAXSI(:,:,:,1),XYZC_FACE_XSIETA(:,:,:,1))

      call node_to_cell_face(QH_XYZN_6_P(:,:,:,2),imax+1,jmax+1,halo_width+1, &
           XYZC_FACE_ETAZTA(:,:,:,2),XYZC_FACE_ZTAXSI(:,:,:,2),XYZC_FACE_XSIETA(:,:,:,2))

      call node_to_cell_face(QH_XYZN_6_P(:,:,:,3),imax+1,jmax+1,halo_width+1, &
           XYZC_FACE_ETAZTA(:,:,:,3),XYZC_FACE_ZTAXSI(:,:,:,3),XYZC_FACE_XSIETA(:,:,:,3))

      call node_to_cell_center(QH_XYZN_6_P(:,:,:,1),imax+1,jmax+1,halo_width+1,XYZC(:,:,:,1))
      call node_to_cell_center(QH_XYZN_6_P(:,:,:,2),imax+1,jmax+1,halo_width+1,XYZC(:,:,:,2))
      call node_to_cell_center(QH_XYZN_6_P(:,:,:,3),imax+1,jmax+1,halo_width+1,XYZC(:,:,:,3))

      call grid_scale(QH_XYZN_6_P(:,:,:,1),QH_XYZN_6_P(:,:,:,2),QH_XYZN_6_P(:,:,:,3), &
                      imax+1,jmax+1,halo_width+1,GRSC)

      do k=BLK%kend_wh+1-halo_width+1,BLK%kend_wh+1
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh+1-halo_width+1)+1
            BLK%XCRD_NODE(i,j,k) = QH_XYZN_6_P(ii,jj,kk+1,1)
            BLK%YCRD_NODE(i,j,k) = QH_XYZN_6_P(ii,jj,kk+1,2)
            BLK%ZCRD_NODE(i,j,k) = QH_XYZN_6_P(ii,jj,kk+1,3)
          enddo
        enddo
      enddo

      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend+1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh-halo_width+1)+1
            BLK%XCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,1)
            BLK%YCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,2)
            BLK%ZCRD_FACE_ETAZTA(i,j,k) = XYZC_FACE_ETAZTA(ii,jj,kk,3)
            BLK%DXSX_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,1)
            BLK%DXSY_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,2)
            BLK%DXSZ_FACE(i,j,k) = XSIX_FACE(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend+1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh-halo_width+1)+1
            BLK%XCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,1)
            BLK%YCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,2)
            BLK%ZCRD_FACE_ZTAXSI(i,j,k) = XYZC_FACE_ZTAXSI(ii,jj,kk,3)
            BLK%DETX_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,1)
            BLK%DETY_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,2)
            BLK%DETZ_FACE(i,j,k) = ETAX_FACE(ii,jj,kk,3)
          enddo
        enddo
      enddo

      do k=BLK%kend_wh+1-halo_width+1,BLK%kend_wh+1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh+1-halo_width+1)+1
            BLK%XCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk+1,1)
            BLK%YCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk+1,2)
            BLK%ZCRD_FACE_XSIETA(i,j,k) = XYZC_FACE_XSIETA(ii,jj,kk+1,3)
            BLK%DZTX_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk+1,1)
            BLK%DZTY_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk+1,2)
            BLK%DZTZ_FACE(i,j,k) = ZTAX_FACE(ii,jj,kk+1,3)
          enddo
        enddo
      enddo

      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh-halo_width+1)+1
            BLK%DXSX(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,1)+XSIX_FACE(ii,jj,kk,1))
            BLK%DXSY(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,2)+XSIX_FACE(ii,jj,kk,2))
            BLK%DXSZ(i,j,k) = 0.5d0*(XSIX_FACE(ii+1,jj,kk,3)+XSIX_FACE(ii,jj,kk,3))
            BLK%DETX(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,1)+ETAX_FACE(ii,jj,kk,1))
            BLK%DETY(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,2)+ETAX_FACE(ii,jj,kk,2))
            BLK%DETZ(i,j,k) = 0.5d0*(ETAX_FACE(ii,jj+1,kk,3)+ETAX_FACE(ii,jj,kk,3))
            BLK%DZTX(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,1)+ZTAX_FACE(ii,jj,kk,1))
            BLK%DZTY(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,2)+ZTAX_FACE(ii,jj,kk,2))
            BLK%DZTZ(i,j,k) = 0.5d0*(ZTAX_FACE(ii,jj,kk+1,3)+ZTAX_FACE(ii,jj,kk,3))
            BLK%XCRD(i,j,k) = XYZC(ii,jj,kk,1)
            BLK%YCRD(i,j,k) = XYZC(ii,jj,kk,2)
            BLK%ZCRD(i,j,k) = XYZC(ii,jj,kk,3)
            BLK%RJCB(i,j,k) = RJCB(ii,jj,kk)
            BLK%GRSC(i,j,k) = GRSC(ii,jj,kk)
          enddo
        enddo
      enddo

    endif

    deallocate(XSIX_FACE)
    deallocate(ETAX_FACE)
    deallocate(ZTAX_FACE)
    deallocate(XYZC)
    deallocate(XYZC_FACE_ETAZTA)
    deallocate(XYZC_FACE_ZTAXSI)
    deallocate(XYZC_FACE_XSIETA)
    deallocate(RJCB)
    deallocate(GRSC)

    deallocate(QH_XYZN_1_M)
    deallocate(QH_XYZN_1_P)
    deallocate(QH_XYZN_2_M)
    deallocate(QH_XYZN_2_P)
    deallocate(QH_XYZN_3_M)
    deallocate(QH_XYZN_3_P)
    deallocate(QH_XYZN_4_M)
    deallocate(QH_XYZN_4_P)
    deallocate(QH_XYZN_5_M)
    deallocate(QH_XYZN_5_P)
    deallocate(QH_XYZN_6_M)
    deallocate(QH_XYZN_6_P)
  end subroutine halo_exchange_metric
!-------------------------------------------------------------------------------
  subroutine halo_exchange_array(BLK,QARR,IFLAG_ARRY,lmax)
    type(block), intent(inout) :: BLK
    real(8)    , intent(inout), dimension(:,:,:,:) :: QARR
    integer    , intent(out)  , dimension(:) :: IFLAG_ARRY
    integer    , intent(in) :: lmax

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: imax,jmax,kmax
    integer :: i,j,k,l,ii,jj,kk

    integer :: iact_xsi,iact_eta,iact_zta

    integer, allocatable, dimension(:) :: IFLAG_TEMP

    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_1_M, QH_TEMP_1_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_2_M, QH_TEMP_2_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_3_M, QH_TEMP_3_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_4_M, QH_TEMP_4_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_5_M, QH_TEMP_5_P
    real(8), allocatable, dimension(:,:,:,:) :: QH_TEMP_6_M, QH_TEMP_6_P

    if (halo_width<=0) return

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainPeriodicFlag_From_DomainID(id_domain,iact_xsi,iact_eta,iact_zta)

    allocate(IFLAG_TEMP(6))

    allocate(QH_TEMP_1_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_1_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_2_M(halo_width,jmax,kmax,lmax))
    allocate(QH_TEMP_2_P(halo_width,jmax,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_3_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_3_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_4_M(imax,halo_width,kmax,lmax))
    allocate(QH_TEMP_4_P(imax,halo_width,kmax,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_5_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_5_P(imax,jmax,halo_width,lmax))
    !-----------------------------------------------
    allocate(QH_TEMP_6_M(imax,jmax,halo_width,lmax))
    allocate(QH_TEMP_6_P(imax,jmax,halo_width,lmax))

    QH_TEMP_1_M(:,:,:,:) = 0.0d0
    QH_TEMP_1_P(:,:,:,:) = 0.0d0
    QH_TEMP_2_M(:,:,:,:) = 0.0d0
    QH_TEMP_2_P(:,:,:,:) = 0.0d0
    QH_TEMP_3_M(:,:,:,:) = 0.0d0
    QH_TEMP_3_P(:,:,:,:) = 0.0d0
    QH_TEMP_4_M(:,:,:,:) = 0.0d0
    QH_TEMP_4_P(:,:,:,:) = 0.0d0
    QH_TEMP_5_M(:,:,:,:) = 0.0d0
    QH_TEMP_5_P(:,:,:,:) = 0.0d0
    QH_TEMP_6_M(:,:,:,:) = 0.0d0
    QH_TEMP_6_P(:,:,:,:) = 0.0d0

    !1
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%ista+halo_width-1
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_1_M(ii,jj,kk,l) = QARR(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    !2
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend-halo_width+1,BLK%iend
            ii = i-(BLK%iend-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_2_M(ii,jj,kk,l) = QARR(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    !3
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jsta+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_3_M(ii,jj,kk,l) = QARR(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    !4
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend-halo_width+1,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend-halo_width+1)+1
            kk = k-BLK%ksta+1
            QH_TEMP_4_M(ii,jj,kk,l) = QARR(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    if (.not.(i2dimension>0)) then
    !5
    do l=1,lmax
      do k=BLK%ksta,BLK%ksta+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QH_TEMP_5_M(ii,jj,kk,l) = QARR(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (.not.(i2dimension>0)) then
    !6
    do l=1,lmax
      do k=BLK%kend-halo_width+1,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend-halo_width+1)+1
            QH_TEMP_6_M(ii,jj,kk,l) = QARR(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
    endif

    IFLAG_TEMP(:) = 0

    !exchange
    call exchange_domain_halo(QH_TEMP_1_M,QH_TEMP_1_P,QH_TEMP_2_M,QH_TEMP_2_P, &
                              QH_TEMP_3_M,QH_TEMP_3_P,QH_TEMP_4_M,QH_TEMP_4_P, &
                              QH_TEMP_5_M,QH_TEMP_5_P,QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                 QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)

    call mpisub_sbsp_exchange_subdomain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                 QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)

    if (.not.(i2dimension>0)) then
      call mpisub_sbsp_exchange_subdomain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                   QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_xsi>0) then
      call mpisub_sbsp_exchange_domain_halo_xsi(QH_TEMP_1_M,QH_TEMP_1_P, &
                                                QH_TEMP_2_M,QH_TEMP_2_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if (iact_eta>0) then
      call mpisub_sbsp_exchange_domain_halo_eta(QH_TEMP_3_M,QH_TEMP_3_P, &
                                                QH_TEMP_4_M,QH_TEMP_4_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    if ((.not.(i2dimension>0)).and.(iact_zta>0)) then
      call mpisub_sbsp_exchange_domain_halo_zta(QH_TEMP_5_M,QH_TEMP_5_P, &
                                                QH_TEMP_6_M,QH_TEMP_6_P, lmax, halo_width, IFLAG_TEMP, 0)
    endif

    IFLAG_ARRY(:) = IFLAG_TEMP(:)

    if (IFLAG_TEMP(1)==1) then
    !1
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            ii = i-BLK%ista_wh+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QARR(i,j,k,l) = QH_TEMP_1_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(2)==1) then
    !2
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            ii = i-(BLK%iend_wh-halo_width+1)+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QARR(i,j,k,l) = QH_TEMP_2_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(3)==1) then
    !3
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta_wh+1
            kk = k-BLK%ksta+1
            QARR(i,j,k,l) = QH_TEMP_3_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if (IFLAG_TEMP(4)==1) then
    !4
    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-(BLK%jend_wh-halo_width+1)+1
            kk = k-BLK%ksta+1
            QARR(i,j,k,l) = QH_TEMP_4_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(5)==1) then
    !5
    do l=1,lmax
      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta_wh+1
            QARR(i,j,k,l) = QH_TEMP_5_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    if ((.not.(i2dimension>0)).and.IFLAG_TEMP(6)==1) then
    !6
    do l=1,lmax
      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-(BLK%kend_wh-halo_width+1)+1
            QARR(i,j,k,l) = QH_TEMP_6_P(ii,jj,kk,l)
          enddo
        enddo
      enddo
    enddo
    endif

    deallocate(IFLAG_TEMP)
    deallocate(QH_TEMP_1_M)
    deallocate(QH_TEMP_1_P)
    deallocate(QH_TEMP_2_M)
    deallocate(QH_TEMP_2_P)
    deallocate(QH_TEMP_3_M)
    deallocate(QH_TEMP_3_P)
    deallocate(QH_TEMP_4_M)
    deallocate(QH_TEMP_4_P)
    deallocate(QH_TEMP_5_M)
    deallocate(QH_TEMP_5_P)
    deallocate(QH_TEMP_6_M)
    deallocate(QH_TEMP_6_P)
  end subroutine halo_exchange_array
!-------------------------------------------------------------------------------
  subroutine halo_init
    integer :: myrank,id_block
    integer :: imax,jmax,kmax

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(IFLAG_PRIM(6))
    allocate(IFLAG_GRAD(6))
    allocate(IFLAG_MTRC(6))
    allocate(IFLAG_RHSS(6))
    allocate(IFLAG_QCS1(6))

    IFLAG_PRIM(:) = 0
    IFLAG_GRAD(:) = 0
    IFLAG_MTRC(:) = 0
    IFLAG_RHSS(:) = 0
    IFLAG_QCS1(:) = 0
  end subroutine halo_init
!-------------------------------------------------------------------------------
  subroutine halo_final
    deallocate(IFLAG_PRIM)
    deallocate(IFLAG_GRAD)
    deallocate(IFLAG_MTRC)
    deallocate(IFLAG_RHSS)
    deallocate(IFLAG_QCS1)
  end subroutine halo_final

end module mdl_halo
