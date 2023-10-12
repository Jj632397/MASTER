module mdl_filter
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_halo
  use mdl_lp_filter
  implicit none

  !input------------------------------------------------------------------------
  integer, private, parameter :: ntimes = 1
  integer, private, parameter :: icompact = 1
  integer, private, parameter, dimension(11) :: odr_xsi = (/    0,     10,     10,     10,     10,     10,     10,     10,     10,     10,    0/)
  integer, private, parameter, dimension(11) :: odr_eta = (/    0,     10,     10,     10,     10,     10,     10,     10,     10,     10,    0/)
  integer, private, parameter, dimension(11) :: odr_zta = (/    0,     10,     10,     10,     10,     10,     10,     10,     10,     10,    0/)
  real(8), private, parameter, dimension(11) :: afs_xsi = (/0.0d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.0d0/)
  real(8), private, parameter, dimension(11) :: afs_eta = (/0.0d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.0d0/)
  real(8), private, parameter, dimension(11) :: afs_zta = (/0.0d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.495d0,0.0d0/)
  !-----------------------------------------------------------------------------

contains

  subroutine filter(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,imax,jmax,kmax
    integer :: n

    integer :: i,j,k,l,ii,jj,kk
    integer :: nhalo1_xsi,nhalo2_xsi
    integer :: nhalo1_eta,nhalo2_eta
    integer :: nhalo1_zta,nhalo2_zta
    integer :: iimax,jjmax
    integer :: ioft,joft,koft

    real(8), allocatable, dimension(:,:) :: QCS1_ALN_XSI
    real(8), allocatable, dimension(:,:) :: QCS1_ALN_ETA
    real(8), allocatable, dimension(:,:) :: QCS1_ALN_ZTA

    integer, allocatable, dimension(:) :: plc_odr_xsi
    integer, allocatable, dimension(:) :: plc_odr_eta
    integer, allocatable, dimension(:) :: plc_odr_zta
    real(8), allocatable, dimension(:) :: plc_afs_xsi
    real(8), allocatable, dimension(:) :: plc_afs_eta
    real(8), allocatable, dimension(:) :: plc_afs_zta

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    do n=1,ntimes

      allocate(plc_odr_xsi(11))
      allocate(plc_odr_eta(11))
      allocate(plc_odr_zta(11))
      allocate(plc_afs_xsi(11))
      allocate(plc_afs_eta(11))
      allocate(plc_afs_zta(11))

      call halo_exchange_qcs(BLK)

      !xsi
      plc_odr_xsi(:) = odr_xsi(:)
      plc_afs_xsi(:) = afs_xsi(:)

      if (IFLAG_QCS1(1)==1) then
        if (halo_width>=2) then
          plc_odr_xsi(2) = 0
          plc_afs_xsi(2) = 0.0d0
        endif
        if (halo_width>=3) then
          plc_odr_xsi(3) = 0
          plc_afs_xsi(3) = 0.0d0
        endif
        if (halo_width>=4) then
          plc_odr_xsi(4) = 0
          plc_afs_xsi(4) = 0.0d0
        endif
        if (halo_width>=5) then
          plc_odr_xsi(5) = 0
          plc_afs_xsi(5) = 0.0d0
        endif
      endif

      if (IFLAG_QCS1(2)==1) then
        if (halo_width>=2) then
          plc_odr_xsi(10) = 0
          plc_afs_xsi(10) = 0.0d0
        endif
        if (halo_width>=3) then
          plc_odr_xsi(9) = 0
          plc_afs_xsi(9) = 0.0d0
        endif
        if (halo_width>=4) then
          plc_odr_xsi(8) = 0
          plc_afs_xsi(8) = 0.0d0
        endif
        if (halo_width>=5) then
          plc_odr_xsi(7) = 0
          plc_afs_xsi(7) = 0.0d0
        endif
      endif

      !eta
      plc_odr_eta(:) = odr_eta(:)
      plc_afs_eta(:) = afs_eta(:)

      if (IFLAG_QCS1(3)==1) then
        if (halo_width>=2) then
          plc_odr_eta(2) = 0
          plc_afs_eta(2) = 0.0d0
        endif
        if (halo_width>=3) then
          plc_odr_eta(3) = 0
          plc_afs_eta(3) = 0.0d0
        endif
        if (halo_width>=4) then
          plc_odr_eta(4) = 0
          plc_afs_eta(4) = 0.0d0
        endif
        if (halo_width>=5) then
          plc_odr_eta(5) = 0
          plc_afs_eta(5) = 0.0d0
        endif
      endif

      if (IFLAG_QCS1(4)==1) then
        if (halo_width>=2) then
          plc_odr_eta(10) = 0
          plc_afs_eta(10) = 0.0d0
        endif
        if (halo_width>=3) then
          plc_odr_eta(9) = 0
          plc_afs_eta(9) = 0.0d0
        endif
        if (halo_width>=4) then
          plc_odr_eta(8) = 0
          plc_afs_eta(8) = 0.0d0
        endif
        if (halo_width>=5) then
          plc_odr_eta(7) = 0
          plc_afs_eta(7) = 0.0d0
        endif
      endif

      !zta
      plc_odr_zta(:) = odr_zta(:)
      plc_afs_zta(:) = afs_zta(:)

      if (IFLAG_QCS1(5)==1) then
        if (halo_width>=2) then
          plc_odr_zta(2) = 0
          plc_afs_zta(2) = 0.0d0
        endif
        if (halo_width>=3) then
          plc_odr_zta(3) = 0
          plc_afs_zta(3) = 0.0d0
        endif
        if (halo_width>=4) then
          plc_odr_zta(4) = 0
          plc_afs_zta(4) = 0.0d0
        endif
        if (halo_width>=5) then
          plc_odr_zta(5) = 0
          plc_afs_zta(5) = 0.0d0
        endif
      endif

      if (IFLAG_QCS1(6)==1) then
        if (halo_width>=2) then
          plc_odr_zta(10) = 0
          plc_afs_zta(10) = 0.0d0
        endif
        if (halo_width>=3) then
          plc_odr_zta(9) = 0
          plc_afs_zta(9) = 0.0d0
        endif
        if (halo_width>=4) then
          plc_odr_zta(8) = 0
          plc_afs_zta(8) = 0.0d0
        endif
        if (halo_width>=5) then
          plc_odr_zta(7) = 0
          plc_afs_zta(7) = 0.0d0
        endif
      endif

      nhalo1_xsi = 0
      nhalo2_xsi = 0
      nhalo1_eta = 0
      nhalo2_eta = 0
      nhalo1_zta = 0
      nhalo2_zta = 0

      if (IFLAG_QCS1(1)==1) nhalo1_xsi = halo_width
      if (IFLAG_QCS1(2)==1) nhalo2_xsi = halo_width
      if (IFLAG_QCS1(3)==1) nhalo1_eta = halo_width
      if (IFLAG_QCS1(4)==1) nhalo2_eta = halo_width
      if (IFLAG_QCS1(5)==1) nhalo1_zta = halo_width
      if (IFLAG_QCS1(6)==1) nhalo2_zta = halo_width

      allocate(QCS1_ALN_XSI(jmax*kmax,imax+nhalo1_xsi+nhalo2_xsi))
      allocate(QCS1_ALN_ETA(imax*kmax,jmax+nhalo1_eta+nhalo2_eta))
      allocate(QCS1_ALN_ZTA(imax*jmax,kmax+nhalo1_zta+nhalo2_zta))


      do l=1,neqns

        !xsi
        ioft = BLK%ista-1
        joft = BLK%jsta-1
        koft = BLK%ksta-1

        if (IFLAG_QCS1(1)==1) ioft = 0

        iimax = jmax*kmax
        jjmax = imax+nhalo1_xsi+nhalo2_xsi

        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + jj
            j = joft + mod(ii-1,jmax)+1
            k = koft + (ii-1)/jmax+1
            QCS1_ALN_XSI(ii,jj) = BLK%QCS1(i,j,k,l)
          enddo
        enddo

        call calc_lp_filter_align(QCS1_ALN_XSI,0,icompact,plc_odr_xsi,plc_afs_xsi,iimax,jjmax)

        do jj=1,imax
          do ii=1,jmax*kmax
            i = BLK%ista-1 + jj
            j = BLK%jsta-1 + mod(ii-1,jmax)+1
            k = BLK%ksta-1 + (ii-1)/jmax+1
            BLK%QCS1(i,j,k,l) = QCS1_ALN_XSI(ii,nhalo1_xsi+jj)
          enddo
        enddo


        !eta
        ioft = BLK%ista-1
        joft = BLK%jsta-1
        koft = BLK%ksta-1

        if (IFLAG_QCS1(3)==1) joft = 0

        iimax = imax*kmax
        jjmax = jmax+nhalo1_eta+nhalo2_eta

        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + jj
            k = koft + (ii-1)/imax+1
            QCS1_ALN_ETA(ii,jj) = BLK%QCS1(i,j,k,l)
          enddo
        enddo

        call calc_lp_filter_align(QCS1_ALN_ETA,0,icompact,plc_odr_eta,plc_afs_eta,iimax,jjmax)

        do jj=1,jmax
          do ii=1,imax*kmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + jj
            k = BLK%ksta-1 + (ii-1)/imax+1
            BLK%QCS1(i,j,k,l) = QCS1_ALN_ETA(ii,nhalo1_eta+jj)
          enddo
        enddo

        !zta
        ioft = BLK%ista-1
        joft = BLK%jsta-1
        koft = BLK%ksta-1

        if (IFLAG_QCS1(5)==1) koft = 0

        iimax = imax*jmax
        jjmax = kmax+nhalo1_zta+nhalo2_zta

        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + (ii-1)/imax+1
            k = koft + jj
            QCS1_ALN_ZTA(ii,jj) = BLK%QCS1(i,j,k,l)
          enddo
        enddo

        call calc_lp_filter_align(QCS1_ALN_ZTA,0,icompact,plc_odr_zta,plc_afs_zta,iimax,jjmax)

        do jj=1,kmax
          do ii=1,imax*jmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + (ii-1)/imax+1
            k = BLK%ksta-1 + jj
            BLK%QCS1(i,j,k,l) = QCS1_ALN_ZTA(ii,nhalo1_zta+jj)
          enddo
        enddo

      enddo

      deallocate(QCS1_ALN_XSI)
      deallocate(QCS1_ALN_ETA)
      deallocate(QCS1_ALN_ZTA)

      deallocate(plc_odr_xsi)
      deallocate(plc_odr_eta)
      deallocate(plc_odr_zta)
      deallocate(plc_afs_xsi)
      deallocate(plc_afs_eta)
      deallocate(plc_afs_zta)
    enddo

  end subroutine filter

end module mdl_filter
