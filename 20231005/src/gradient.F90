module mdl_grad
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_fdm1
  use mdl_halo
  implicit none

  integer, private, parameter :: icompact = 0
  integer, private, parameter, dimension(5) :: isc_xsi = (/1,1,1,1,1/)
  integer, private, parameter, dimension(5) :: isc_eta = (/1,1,1,1,1/)
  integer, private, parameter, dimension(5) :: isc_zta = (/1,1,1,1,1/)

contains

  subroutine gradient(BLK)
    type(block), intent(inout) :: BLK

    integer :: i,j,k,l
    integer :: imax,jmax,kmax,lmax
    integer :: imax_wh,jmax_wh,kmax_wh
    integer :: nhalo1_xsi,nhalo2_xsi
    integer :: nhalo1_eta,nhalo2_eta
    integer :: nhalo1_zta,nhalo2_zta
    integer :: ioft,joft,koft
    integer :: ii,jj,kk
    integer :: iimax,jjmax

    real(8), allocatable, dimension(:,:,:) :: QQ_ALN

    real(8), allocatable, dimension(:,:,:,:) :: DQXS,DQET,DQZT

    real(8), allocatable, dimension(:,:,:) :: RRJCB

    imax = BLK%imax
    jmax = BLK%jmax
    kmax = BLK%kmax

    imax_wh = BLK%imax_wh
    jmax_wh = BLK%jmax_wh
    kmax_wh = BLK%kmax_wh

    lmax = neqns+1

    nhalo1_xsi = 0
    nhalo2_xsi = 0
    nhalo1_eta = 0
    nhalo2_eta = 0
    nhalo1_zta = 0
    nhalo2_zta = 0

    if (IFLAG_PRIM(1)==1) nhalo1_xsi = halo_width
    if (IFLAG_PRIM(2)==1) nhalo2_xsi = halo_width
    if (IFLAG_PRIM(3)==1) nhalo1_eta = halo_width
    if (IFLAG_PRIM(4)==1) nhalo2_eta = halo_width
    if (IFLAG_PRIM(5)==1) nhalo1_zta = halo_width
    if (IFLAG_PRIM(6)==1) nhalo2_zta = halo_width

    allocate(DQXS(imax_wh,jmax_wh,kmax_wh,lmax))
    allocate(DQET(imax_wh,jmax_wh,kmax_wh,lmax))
    allocate(DQZT(imax_wh,jmax_wh,kmax_wh,lmax))
    DQXS(:,:,:,:) = 0.0d0 !Never forget to initialize
    DQET(:,:,:,:) = 0.0d0 !
    DQZT(:,:,:,:) = 0.0d0 !

    allocate(RRJCB(imax_wh,jmax_wh,kmax_wh))
    RRJCB(:,:,:) = 1.0d0 !Never forget to initialize by 1

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          RRJCB(i,j,k) = 1.0d0/BLK%RJCB(i,j,k)
        enddo
      enddo
    enddo

    !Xsi derivative
    iimax = jmax*kmax
    jjmax = imax+nhalo1_xsi+nhalo2_xsi

    allocate(QQ_ALN(iimax,jjmax,lmax))

    ioft = BLK%ista-1
    joft = BLK%jsta-1
    koft = BLK%ksta-1

    if (IFLAG_PRIM(1)==1) ioft = 0

    do jj=1,jjmax
      do ii=1,iimax
        i = ioft + jj
        j = joft + mod(ii-1,jmax)+1
        k = koft + (ii-1)/jmax+1
        QQ_ALN(ii,jj,1)       = BLK%TMPR(i,j,k)
        QQ_ALN(ii,jj,2)       = BLK%ULCT(i,j,k)
        QQ_ALN(ii,jj,3)       = BLK%VLCT(i,j,k)
        QQ_ALN(ii,jj,4)       = BLK%WLCT(i,j,k)
        QQ_ALN(ii,jj,5)       = BLK%PRSS(i,j,k)
        QQ_ALN(ii,jj,neqns+1) = BLK%DNST(i,j,k)
      enddo
    enddo

    do l=1,neqns-5
      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + jj
          j = joft + mod(ii-1,jmax)+1
          k = koft + (ii-1)/jmax+1
          QQ_ALN(ii,jj,5+l) = BLK%YSPC(i,j,k,l)
        enddo
      enddo
    enddo

    do l=1,lmax
      call calc_fdm_1st_align(QQ_ALN(:,:,l),QQ_ALN(:,:,l),0,icompact,isc_xsi,iimax,jjmax)
    enddo

    do l=1,lmax
      do jj=1,imax
        do ii=1,jmax*kmax
          i = BLK%ista-1 + jj
          j = BLK%jsta-1 + mod(ii-1,jmax)+1
          k = BLK%ksta-1 + (ii-1)/jmax+1
          DQXS(i,j,k,l) = QQ_ALN(ii,jj+nhalo1_xsi,l)
        enddo
      enddo
    enddo

    deallocate(QQ_ALN)

    !Eta derivative
    iimax = imax*kmax
    jjmax = jmax+nhalo1_eta+nhalo2_eta

    allocate(QQ_ALN(iimax,jjmax,lmax))

    ioft = BLK%ista-1
    joft = BLK%jsta-1
    koft = BLK%ksta-1

    if (IFLAG_PRIM(3)==1) joft = 0

    do jj=1,jjmax
      do ii=1,iimax
        i = ioft + mod(ii-1,imax)+1
        j = joft + jj
        k = koft + (ii-1)/imax+1
        QQ_ALN(ii,jj,1)       = BLK%TMPR(i,j,k)
        QQ_ALN(ii,jj,2)       = BLK%ULCT(i,j,k)
        QQ_ALN(ii,jj,3)       = BLK%VLCT(i,j,k)
        QQ_ALN(ii,jj,4)       = BLK%WLCT(i,j,k)
        QQ_ALN(ii,jj,5)       = BLK%PRSS(i,j,k)
        QQ_ALN(ii,jj,neqns+1) = BLK%DNST(i,j,k)
      enddo
    enddo

    do l=1,neqns-5
      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + jj
          k = koft + (ii-1)/imax+1
          QQ_ALN(ii,jj,5+l) = BLK%YSPC(i,j,k,l)
        enddo
      enddo
    enddo

    do l=1,lmax
      call calc_fdm_1st_align(QQ_ALN(:,:,l),QQ_ALN(:,:,l),0,icompact,isc_eta,iimax,jjmax)
    enddo

    do l=1,lmax
      do jj=1,jmax
        do ii=1,imax*kmax
          i = BLK%ista-1 + mod(ii-1,imax)+1
          j = BLK%jsta-1 + jj
          k = BLK%ksta-1 + (ii-1)/imax+1
          DQET(i,j,k,l) = QQ_ALN(ii,jj+nhalo1_eta,l)
        enddo
      enddo
    enddo

    deallocate(QQ_ALN)

    !Zta derivative
    if (.not.(i2dimension>0)) then

    iimax = imax*jmax
    jjmax = kmax+nhalo1_zta+nhalo2_zta

    allocate(QQ_ALN(iimax,jjmax,lmax))

    ioft = BLK%ista-1
    joft = BLK%jsta-1
    koft = BLK%ksta-1

    if (IFLAG_PRIM(5)==1) koft = 0

    do jj=1,jjmax
      do ii=1,iimax
        i = ioft + mod(ii-1,imax)+1
        j = joft + (ii-1)/imax+1
        k = koft + jj
        QQ_ALN(ii,jj,1)       = BLK%TMPR(i,j,k)
        QQ_ALN(ii,jj,2)       = BLK%ULCT(i,j,k)
        QQ_ALN(ii,jj,3)       = BLK%VLCT(i,j,k)
        QQ_ALN(ii,jj,4)       = BLK%WLCT(i,j,k)
        QQ_ALN(ii,jj,5)       = BLK%PRSS(i,j,k)
        QQ_ALN(ii,jj,neqns+1) = BLK%DNST(i,j,k)
      enddo
    enddo

    do l=1,neqns-5
      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + (ii-1)/imax+1
          k = koft + jj
          QQ_ALN(ii,jj,5+l) = BLK%YSPC(i,j,k,l)
        enddo
      enddo
    enddo

    do l=1,lmax
      call calc_fdm_1st_align(QQ_ALN(:,:,l),QQ_ALN(:,:,l),0,icompact,isc_zta,iimax,jjmax)
    enddo

    do l=1,lmax
      do jj=1,kmax
        do ii=1,imax*jmax
          i = BLK%ista-1 + mod(ii-1,imax)+1
          j = BLK%jsta-1 + (ii-1)/imax+1
          k = BLK%ksta-1 + jj
          DQZT(i,j,k,l) = QQ_ALN(ii,jj+nhalo1_zta,l)
        enddo
      enddo
    enddo

    deallocate(QQ_ALN)

    else

      DQZT(:,:,:,:) = 0.0d0

    endif

    !--------------------------------------------------------------------------------------
    BLK%DTMX = (BLK%DXSX*DQXS(:,:,:,1)+BLK%DETX*DQET(:,:,:,1)+BLK%DZTX*DQZT(:,:,:,1))*RRJCB
    BLK%DTMY = (BLK%DXSY*DQXS(:,:,:,1)+BLK%DETY*DQET(:,:,:,1)+BLK%DZTY*DQZT(:,:,:,1))*RRJCB
    BLK%DTMZ = (BLK%DXSZ*DQXS(:,:,:,1)+BLK%DETZ*DQET(:,:,:,1)+BLK%DZTZ*DQZT(:,:,:,1))*RRJCB
    !--------------------------------------------------------------------------------------
    BLK%DULX = (BLK%DXSX*DQXS(:,:,:,2)+BLK%DETX*DQET(:,:,:,2)+BLK%DZTX*DQZT(:,:,:,2))*RRJCB
    BLK%DULY = (BLK%DXSY*DQXS(:,:,:,2)+BLK%DETY*DQET(:,:,:,2)+BLK%DZTY*DQZT(:,:,:,2))*RRJCB
    BLK%DULZ = (BLK%DXSZ*DQXS(:,:,:,2)+BLK%DETZ*DQET(:,:,:,2)+BLK%DZTZ*DQZT(:,:,:,2))*RRJCB
    !--------------------------------------------------------------------------------------
    BLK%DVLX = (BLK%DXSX*DQXS(:,:,:,3)+BLK%DETX*DQET(:,:,:,3)+BLK%DZTX*DQZT(:,:,:,3))*RRJCB
    BLK%DVLY = (BLK%DXSY*DQXS(:,:,:,3)+BLK%DETY*DQET(:,:,:,3)+BLK%DZTY*DQZT(:,:,:,3))*RRJCB
    BLK%DVLZ = (BLK%DXSZ*DQXS(:,:,:,3)+BLK%DETZ*DQET(:,:,:,3)+BLK%DZTZ*DQZT(:,:,:,3))*RRJCB
    !--------------------------------------------------------------------------------------
    BLK%DWLX = (BLK%DXSX*DQXS(:,:,:,4)+BLK%DETX*DQET(:,:,:,4)+BLK%DZTX*DQZT(:,:,:,4))*RRJCB
    BLK%DWLY = (BLK%DXSY*DQXS(:,:,:,4)+BLK%DETY*DQET(:,:,:,4)+BLK%DZTY*DQZT(:,:,:,4))*RRJCB
    BLK%DWLZ = (BLK%DXSZ*DQXS(:,:,:,4)+BLK%DETZ*DQET(:,:,:,4)+BLK%DZTZ*DQZT(:,:,:,4))*RRJCB
    !--------------------------------------------------------------------------------------
    BLK%DPRX = (BLK%DXSX*DQXS(:,:,:,5)+BLK%DETX*DQET(:,:,:,5)+BLK%DZTX*DQZT(:,:,:,5))*RRJCB
    BLK%DPRY = (BLK%DXSY*DQXS(:,:,:,5)+BLK%DETY*DQET(:,:,:,5)+BLK%DZTY*DQZT(:,:,:,5))*RRJCB
    BLK%DPRZ = (BLK%DXSZ*DQXS(:,:,:,5)+BLK%DETZ*DQET(:,:,:,5)+BLK%DZTZ*DQZT(:,:,:,5))*RRJCB
    !--------------------------------------------------------------------------------------
    BLK%DDNX = (BLK%DXSX*DQXS(:,:,:,neqns+1)+BLK%DETX*DQET(:,:,:,neqns+1)+BLK%DZTX*DQZT(:,:,:,neqns+1))*RRJCB
    BLK%DDNY = (BLK%DXSY*DQXS(:,:,:,neqns+1)+BLK%DETY*DQET(:,:,:,neqns+1)+BLK%DZTY*DQZT(:,:,:,neqns+1))*RRJCB
    BLK%DDNZ = (BLK%DXSZ*DQXS(:,:,:,neqns+1)+BLK%DETZ*DQET(:,:,:,neqns+1)+BLK%DZTZ*DQZT(:,:,:,neqns+1))*RRJCB

    do l=1,neqns-5
      BLK%DYSX(:,:,:,l) = (BLK%DXSX*DQXS(:,:,:,5+l)+BLK%DETX*DQET(:,:,:,5+l)+BLK%DZTX*DQZT(:,:,:,5+l))*RRJCB
      BLK%DYSY(:,:,:,l) = (BLK%DXSY*DQXS(:,:,:,5+l)+BLK%DETY*DQET(:,:,:,5+l)+BLK%DZTY*DQZT(:,:,:,5+l))*RRJCB
      BLK%DYSZ(:,:,:,l) = (BLK%DXSZ*DQXS(:,:,:,5+l)+BLK%DETZ*DQET(:,:,:,5+l)+BLK%DZTZ*DQZT(:,:,:,5+l))*RRJCB
    enddo

    deallocate(DQXS)
    deallocate(DQET)
    deallocate(DQZT)
    deallocate(RRJCB)
  end subroutine gradient

end module mdl_grad
