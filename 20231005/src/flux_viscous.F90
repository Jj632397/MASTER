module mdl_flux_viscous
  use mdl_mpisub_sbsp
  use mdl_decompo
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_fdm1
  use mdl_halo
  use mdl_tbl_sa
  implicit none

  integer, private, parameter :: icompact = 0
  integer, private, parameter, dimension(5) :: isc_xsi = (/1,1,1,1,1/)
  integer, private, parameter, dimension(5) :: isc_eta = (/1,1,1,1,1/)
  integer, private, parameter, dimension(5) :: isc_zta = (/1,1,1,1,1/)
  !-----------------------------------------------------------------------------
  integer, private, parameter :: neqn_passive_scalar_1 = 6  ! Set the index of scalar transport equation.
  integer, private, parameter :: neqn_passive_scalar_2 = 7  ! Then, scolar diffusion is activated.
  integer, private, parameter :: neqn_passive_scalar_3 = 0  ! If diffusion isn't needed, set 0.

contains

  subroutine flux_viscous(BLK,RHS)
    type(block), intent(inout) :: BLK
    real(8), intent(inout), dimension(:,:,:,:) :: RHS

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: iimax,jjmax

    real(8), allocatable, dimension(:,:)   :: DNST_ALN
    real(8), allocatable, dimension(:,:)   :: ULCT_ALN
    real(8), allocatable, dimension(:,:)   :: VLCT_ALN
    real(8), allocatable, dimension(:,:)   :: WLCT_ALN
    real(8), allocatable, dimension(:,:,:) :: YSPC_ALN

    real(8), allocatable, dimension(:,:)   :: DULX_ALN, DULY_ALN, DULZ_ALN
    real(8), allocatable, dimension(:,:)   :: DVLX_ALN, DVLY_ALN, DVLZ_ALN
    real(8), allocatable, dimension(:,:)   :: DWLX_ALN, DWLY_ALN, DWLZ_ALN
    real(8), allocatable, dimension(:,:)   :: DTMX_ALN, DTMY_ALN, DTMZ_ALN
    real(8), allocatable, dimension(:,:,:) :: DYSX_ALN, DYSY_ALN, DYSZ_ALN

    real(8), allocatable, dimension(:,:)   :: VSSV_ALN, VSST_ALN
    real(8), allocatable, dimension(:,:)   :: VSBV_ALN, VSBT_ALN
    real(8), allocatable, dimension(:,:)   :: CKPV_ALN, CKPT_ALN

    real(8), allocatable, dimension(:,:)   :: DXSX_ALN
    real(8), allocatable, dimension(:,:)   :: DXSY_ALN
    real(8), allocatable, dimension(:,:)   :: DXSZ_ALN

    real(8), allocatable, dimension(:,:,:) :: FLUX_ALN
    
    real(8), allocatable, dimension(:,:)   :: TMPR_ALN  !! FURUSAWA

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)


    call calc_xsi
    call calc_eta
    if (.not.(i2dimension>0)) then
      call calc_zta
    endif

  contains
!-------------------------------------------------------------------------------
    subroutine calc_xsi
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft

      nhalo1 = 0
      nhalo2 = 0

      if (IFLAG_PRIM(1)==1.and.IFLAG_GRAD(1)==1) nhalo1 = halo_width
      if (IFLAG_PRIM(2)==1.and.IFLAG_GRAD(2)==1) nhalo2 = halo_width

      iimax = jmax*kmax
      jjmax = imax+nhalo1+nhalo2

      ioft = BLK%ista-1
      joft = BLK%jsta-1
      koft = BLK%ksta-1

      if (IFLAG_PRIM(1)==1.and.IFLAG_GRAD(1)==1) ioft = 0

      allocate(DNST_ALN(iimax,jjmax))
      allocate(ULCT_ALN(iimax,jjmax))
      allocate(VLCT_ALN(iimax,jjmax))
      allocate(WLCT_ALN(iimax,jjmax))
      if (neqns-5>0) then
        allocate(YSPC_ALN(iimax,jjmax,neqns-5))
      endif

      allocate(DULX_ALN(iimax,jjmax))
      allocate(DULY_ALN(iimax,jjmax))
      allocate(DULZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DVLX_ALN(iimax,jjmax))
      allocate(DVLY_ALN(iimax,jjmax))
      allocate(DVLZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DWLX_ALN(iimax,jjmax))
      allocate(DWLY_ALN(iimax,jjmax))
      allocate(DWLZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DTMX_ALN(iimax,jjmax))
      allocate(DTMY_ALN(iimax,jjmax))
      allocate(DTMZ_ALN(iimax,jjmax))
      !------------------------------
      if (neqns-5>0) then
        allocate(DYSX_ALN(iimax,jjmax,neqns-5))
        allocate(DYSY_ALN(iimax,jjmax,neqns-5))
        allocate(DYSZ_ALN(iimax,jjmax,neqns-5))
      endif
      !------------------------------
      allocate(VSSV_ALN(iimax,jjmax))
      allocate(VSST_ALN(iimax,jjmax))
      allocate(VSBV_ALN(iimax,jjmax))
      allocate(VSBT_ALN(iimax,jjmax))
      allocate(CKPV_ALN(iimax,jjmax))
      allocate(CKPT_ALN(iimax,jjmax))
      !------------------------------
      allocate(DXSX_ALN(iimax,jjmax))
      allocate(DXSY_ALN(iimax,jjmax))
      allocate(DXSZ_ALN(iimax,jjmax))

      allocate(FLUX_ALN(iimax,jjmax,neqns))
      
      allocate(TMPR_ALN(iimax,jjmax))   !! FURUSAWA

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + jj
          j = joft + mod(ii-1,jmax)+1
          k = koft + (ii-1)/jmax+1

          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          ULCT_ALN(ii,jj) = BLK%ULCT(i,j,k)
          VLCT_ALN(ii,jj) = BLK%VLCT(i,j,k)
          WLCT_ALN(ii,jj) = BLK%WLCT(i,j,k)
          !--------------------------------
          DULX_ALN(ii,jj) = BLK%DULX(i,j,k)
          DULY_ALN(ii,jj) = BLK%DULY(i,j,k)
          DULZ_ALN(ii,jj) = BLK%DULZ(i,j,k)
          !--------------------------------
          DVLX_ALN(ii,jj) = BLK%DVLX(i,j,k)
          DVLY_ALN(ii,jj) = BLK%DVLY(i,j,k)
          DVLZ_ALN(ii,jj) = BLK%DVLZ(i,j,k)
          !--------------------------------
          DWLX_ALN(ii,jj) = BLK%DWLX(i,j,k)
          DWLY_ALN(ii,jj) = BLK%DWLY(i,j,k)
          DWLZ_ALN(ii,jj) = BLK%DWLZ(i,j,k)
          !--------------------------------
          DTMX_ALN(ii,jj) = BLK%DTMX(i,j,k)
          DTMY_ALN(ii,jj) = BLK%DTMY(i,j,k)
          DTMZ_ALN(ii,jj) = BLK%DTMZ(i,j,k)
          !--------------------------------
          VSSV_ALN(ii,jj) = BLK%VSSV(i,j,k)
          VSST_ALN(ii,jj) = BLK%VSST(i,j,k)
          VSBV_ALN(ii,jj) = BLK%VSBV(i,j,k)
          VSBT_ALN(ii,jj) = BLK%VSBT(i,j,k)
          CKPV_ALN(ii,jj) = BLK%CKPV(i,j,k)
          CKPT_ALN(ii,jj) = BLK%CKPT(i,j,k)
          !--------------------------------
          DXSX_ALN(ii,jj) = BLK%DXSX(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DXSY(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DXSZ(i,j,k)
          !--------------------------------
          TMPR_ALN(ii,jj) = BLK%TMPR(i,j,k)  !! FURUSAWA

        enddo
      enddo

      do l=1,neqns-5
        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + jj
            j = joft + mod(ii-1,jmax)+1
            k = koft + (ii-1)/jmax+1
            YSPC_ALN(ii,jj,l) = BLK%YSPC(i,j,k,l)
            !------------------------------------
            DYSX_ALN(ii,jj,l) = BLK%DYSX(i,j,k,l)
            DYSY_ALN(ii,jj,l) = BLK%DYSY(i,j,k,l)
            DYSZ_ALN(ii,jj,l) = BLK%DYSZ(i,j,k,l)
          enddo
        enddo
      enddo


      call calc_viscous


      do l=1,neqns
        if (    (l==2).or.(l==3).or.(l==4).or.(l==5) &
            .or.(iact_tbl_sa>0.and.l==neqn_tbl_sa)   &
            .or.(l==neqn_passive_scalar_1)           &
            .or.(l==neqn_passive_scalar_2)           &
            .or.(l==neqn_passive_scalar_3)            ) then

          call calc_fdm_1st_align(FLUX_ALN(:,:,l),FLUX_ALN(:,:,l),0,icompact,isc_xsi,iimax,jjmax)

        endif
      enddo


      do l=1,neqns
        if (    (l==2).or.(l==3).or.(l==4).or.(l==5) &
            .or.(iact_tbl_sa>0.and.l==neqn_tbl_sa)   &
            .or.(l==neqn_passive_scalar_1)           &
            .or.(l==neqn_passive_scalar_2)           &
            .or.(l==neqn_passive_scalar_3)            ) then

          do jj=1,imax
            do ii=1,jmax*kmax
              i = BLK%ista-1 + jj
              j = BLK%jsta-1 + mod(ii-1,jmax)+1
              k = BLK%ksta-1 + (ii-1)/jmax+1
              RHS(i,j,k,l) = RHS(i,j,k,l) - FLUX_ALN(ii,nhalo1+jj,l)
            enddo
          enddo

        endif
      enddo

      deallocate(DNST_ALN)
      deallocate(ULCT_ALN)
      deallocate(VLCT_ALN)
      deallocate(WLCT_ALN)
      if (neqns-5>0) then
        deallocate(YSPC_ALN)
      endif
      !-------------------
      deallocate(DULX_ALN)
      deallocate(DULY_ALN)
      deallocate(DULZ_ALN)
      !-------------------
      deallocate(DVLX_ALN)
      deallocate(DVLY_ALN)
      deallocate(DVLZ_ALN)
      !-------------------
      deallocate(DWLX_ALN)
      deallocate(DWLY_ALN)
      deallocate(DWLZ_ALN)
      !-------------------
      deallocate(DTMX_ALN)
      deallocate(DTMY_ALN)
      deallocate(DTMZ_ALN)
      !-------------------
      if (neqns-5>0) then
        deallocate(DYSX_ALN)
        deallocate(DYSY_ALN)
        deallocate(DYSZ_ALN)
      endif
      !-------------------
      deallocate(VSSV_ALN)
      deallocate(VSST_ALN)
      deallocate(VSBV_ALN)
      deallocate(VSBT_ALN)
      deallocate(CKPV_ALN)
      deallocate(CKPT_ALN)
      !-------------------
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      !-------------------
      deallocate(FLUX_ALN)
      !-------------------
      deallocate(TMPR_ALN)  !! FURUSAWA
    end subroutine calc_xsi
!-------------------------------------------------------------------------------
    subroutine calc_eta
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft

      nhalo1 = 0
      nhalo2 = 0

      if (IFLAG_PRIM(3)==1.and.IFLAG_GRAD(3)==1) nhalo1 = halo_width
      if (IFLAG_PRIM(4)==1.and.IFLAG_GRAD(4)==1) nhalo2 = halo_width

      iimax = imax*kmax
      jjmax = jmax+nhalo1+nhalo2

      ioft = BLK%ista-1
      joft = BLK%jsta-1
      koft = BLK%ksta-1

      if (IFLAG_PRIM(3)==1.and.IFLAG_GRAD(3)==1) joft = 0

      allocate(DNST_ALN(iimax,jjmax))
      allocate(ULCT_ALN(iimax,jjmax))
      allocate(VLCT_ALN(iimax,jjmax))
      allocate(WLCT_ALN(iimax,jjmax))
      if (neqns-5>0) then
        allocate(YSPC_ALN(iimax,jjmax,neqns-5))
      endif

      allocate(DULX_ALN(iimax,jjmax))
      allocate(DULY_ALN(iimax,jjmax))
      allocate(DULZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DVLX_ALN(iimax,jjmax))
      allocate(DVLY_ALN(iimax,jjmax))
      allocate(DVLZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DWLX_ALN(iimax,jjmax))
      allocate(DWLY_ALN(iimax,jjmax))
      allocate(DWLZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DTMX_ALN(iimax,jjmax))
      allocate(DTMY_ALN(iimax,jjmax))
      allocate(DTMZ_ALN(iimax,jjmax))
      !------------------------------
      if (neqns-5>0) then
        allocate(DYSX_ALN(iimax,jjmax,neqns-5))
        allocate(DYSY_ALN(iimax,jjmax,neqns-5))
        allocate(DYSZ_ALN(iimax,jjmax,neqns-5))
      endif
      !------------------------------
      allocate(VSSV_ALN(iimax,jjmax))
      allocate(VSST_ALN(iimax,jjmax))
      allocate(VSBV_ALN(iimax,jjmax))
      allocate(VSBT_ALN(iimax,jjmax))
      allocate(CKPV_ALN(iimax,jjmax))
      allocate(CKPT_ALN(iimax,jjmax))
      !------------------------------
      allocate(DXSX_ALN(iimax,jjmax))
      allocate(DXSY_ALN(iimax,jjmax))
      allocate(DXSZ_ALN(iimax,jjmax))

      allocate(FLUX_ALN(iimax,jjmax,neqns))

      allocate(TMPR_ALN(iimax,jjmax))   !! FURUSAWA
      
      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + jj
          k = koft + (ii-1)/imax+1

          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          ULCT_ALN(ii,jj) = BLK%ULCT(i,j,k)
          VLCT_ALN(ii,jj) = BLK%VLCT(i,j,k)
          WLCT_ALN(ii,jj) = BLK%WLCT(i,j,k)
          !--------------------------------
          DULX_ALN(ii,jj) = BLK%DULX(i,j,k)
          DULY_ALN(ii,jj) = BLK%DULY(i,j,k)
          DULZ_ALN(ii,jj) = BLK%DULZ(i,j,k)
          !--------------------------------
          DVLX_ALN(ii,jj) = BLK%DVLX(i,j,k)
          DVLY_ALN(ii,jj) = BLK%DVLY(i,j,k)
          DVLZ_ALN(ii,jj) = BLK%DVLZ(i,j,k)
          !--------------------------------
          DWLX_ALN(ii,jj) = BLK%DWLX(i,j,k)
          DWLY_ALN(ii,jj) = BLK%DWLY(i,j,k)
          DWLZ_ALN(ii,jj) = BLK%DWLZ(i,j,k)
          !--------------------------------
          DTMX_ALN(ii,jj) = BLK%DTMX(i,j,k)
          DTMY_ALN(ii,jj) = BLK%DTMY(i,j,k)
          DTMZ_ALN(ii,jj) = BLK%DTMZ(i,j,k)
          !--------------------------------
          VSSV_ALN(ii,jj) = BLK%VSSV(i,j,k)
          VSST_ALN(ii,jj) = BLK%VSST(i,j,k)
          VSBV_ALN(ii,jj) = BLK%VSBV(i,j,k)
          VSBT_ALN(ii,jj) = BLK%VSBT(i,j,k)
          CKPV_ALN(ii,jj) = BLK%CKPV(i,j,k)
          CKPT_ALN(ii,jj) = BLK%CKPT(i,j,k)
          !--------------------------------
          DXSX_ALN(ii,jj) = BLK%DETX(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DETY(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DETZ(i,j,k)
          !--------------------------------
          TMPR_ALN(ii,jj) = BLK%TMPR(i,j,k)  !! FURUSAWA
        enddo
      enddo

      do l=1,neqns-5
        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + jj
            k = koft + (ii-1)/imax+1
            YSPC_ALN(ii,jj,l) = BLK%YSPC(i,j,k,l)
            !------------------------------------
            DYSX_ALN(ii,jj,l) = BLK%DYSX(i,j,k,l)
            DYSY_ALN(ii,jj,l) = BLK%DYSY(i,j,k,l)
            DYSZ_ALN(ii,jj,l) = BLK%DYSZ(i,j,k,l)
          enddo
        enddo
      enddo


      call calc_viscous


      do l=1,neqns
        if (    (l==2).or.(l==3).or.(l==4).or.(l==5) &
            .or.(iact_tbl_sa>0.and.l==neqn_tbl_sa)   &
            .or.(l==neqn_passive_scalar_1)           &
            .or.(l==neqn_passive_scalar_2)           &
            .or.(l==neqn_passive_scalar_3)            ) then

          call calc_fdm_1st_align(FLUX_ALN(:,:,l),FLUX_ALN(:,:,l),0,icompact,isc_eta,iimax,jjmax)

        endif
      enddo


      do l=1,neqns
        if (    (l==2).or.(l==3).or.(l==4).or.(l==5) &
            .or.(iact_tbl_sa>0.and.l==neqn_tbl_sa)   &
            .or.(l==neqn_passive_scalar_1)           &
            .or.(l==neqn_passive_scalar_2)           &
            .or.(l==neqn_passive_scalar_3)            ) then

          do jj=1,jmax
            do ii=1,imax*kmax
              i = BLK%ista-1 + mod(ii-1,imax)+1
              j = BLK%jsta-1 + jj
              k = BLK%ksta-1 + (ii-1)/imax+1
              RHS(i,j,k,l) = RHS(i,j,k,l) - FLUX_ALN(ii,nhalo1+jj,l)
            enddo
          enddo

        endif
      enddo

      deallocate(DNST_ALN)
      deallocate(ULCT_ALN)
      deallocate(VLCT_ALN)
      deallocate(WLCT_ALN)
      if (neqns-5>0) then
        deallocate(YSPC_ALN)
      endif
      !-------------------
      deallocate(DULX_ALN)
      deallocate(DULY_ALN)
      deallocate(DULZ_ALN)
      !-------------------
      deallocate(DVLX_ALN)
      deallocate(DVLY_ALN)
      deallocate(DVLZ_ALN)
      !-------------------
      deallocate(DWLX_ALN)
      deallocate(DWLY_ALN)
      deallocate(DWLZ_ALN)
      !-------------------
      deallocate(DTMX_ALN)
      deallocate(DTMY_ALN)
      deallocate(DTMZ_ALN)
      !-------------------
      if (neqns-5>0) then
        deallocate(DYSX_ALN)
        deallocate(DYSY_ALN)
        deallocate(DYSZ_ALN)
      endif
      !-------------------
      deallocate(VSSV_ALN)
      deallocate(VSST_ALN)
      deallocate(VSBV_ALN)
      deallocate(VSBT_ALN)
      deallocate(CKPV_ALN)
      deallocate(CKPT_ALN)
      !-------------------
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      !-------------------
      deallocate(FLUX_ALN)
      !-------------------
      deallocate(TMPR_ALN)
    end subroutine calc_eta
!-------------------------------------------------------------------------------
    subroutine calc_zta
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft

      nhalo1 = 0
      nhalo2 = 0

      if (IFLAG_PRIM(5)==1.and.IFLAG_GRAD(5)==1) nhalo1 = halo_width
      if (IFLAG_PRIM(6)==1.and.IFLAG_GRAD(6)==1) nhalo2 = halo_width

      iimax = imax*jmax
      jjmax = kmax+nhalo1+nhalo2

      ioft = BLK%ista-1
      joft = BLK%jsta-1
      koft = BLK%ksta-1

      if (IFLAG_PRIM(5)==1.and.IFLAG_GRAD(5)==1) koft = 0

      allocate(DNST_ALN(iimax,jjmax))
      allocate(ULCT_ALN(iimax,jjmax))
      allocate(VLCT_ALN(iimax,jjmax))
      allocate(WLCT_ALN(iimax,jjmax))
      if (neqns-5>0) then
        allocate(YSPC_ALN(iimax,jjmax,neqns-5))
      endif

      allocate(DULX_ALN(iimax,jjmax))
      allocate(DULY_ALN(iimax,jjmax))
      allocate(DULZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DVLX_ALN(iimax,jjmax))
      allocate(DVLY_ALN(iimax,jjmax))
      allocate(DVLZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DWLX_ALN(iimax,jjmax))
      allocate(DWLY_ALN(iimax,jjmax))
      allocate(DWLZ_ALN(iimax,jjmax))
      !------------------------------
      allocate(DTMX_ALN(iimax,jjmax))
      allocate(DTMY_ALN(iimax,jjmax))
      allocate(DTMZ_ALN(iimax,jjmax))
      !------------------------------
      if (neqns-5>0) then
        allocate(DYSX_ALN(iimax,jjmax,neqns-5))
        allocate(DYSY_ALN(iimax,jjmax,neqns-5))
        allocate(DYSZ_ALN(iimax,jjmax,neqns-5))
      endif
      !------------------------------
      allocate(VSSV_ALN(iimax,jjmax))
      allocate(VSST_ALN(iimax,jjmax))
      allocate(VSBV_ALN(iimax,jjmax))
      allocate(VSBT_ALN(iimax,jjmax))
      allocate(CKPV_ALN(iimax,jjmax))
      allocate(CKPT_ALN(iimax,jjmax))
      !------------------------------
      allocate(DXSX_ALN(iimax,jjmax))
      allocate(DXSY_ALN(iimax,jjmax))
      allocate(DXSZ_ALN(iimax,jjmax))

      allocate(FLUX_ALN(iimax,jjmax,neqns))

      allocate(TMPR_ALN(iimax,jjmax))   !! FURUSAWA

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + (ii-1)/imax+1
          k = koft + jj

          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          ULCT_ALN(ii,jj) = BLK%ULCT(i,j,k)
          VLCT_ALN(ii,jj) = BLK%VLCT(i,j,k)
          WLCT_ALN(ii,jj) = BLK%WLCT(i,j,k)
          !--------------------------------
          DULX_ALN(ii,jj) = BLK%DULX(i,j,k)
          DULY_ALN(ii,jj) = BLK%DULY(i,j,k)
          DULZ_ALN(ii,jj) = BLK%DULZ(i,j,k)
          !--------------------------------
          DVLX_ALN(ii,jj) = BLK%DVLX(i,j,k)
          DVLY_ALN(ii,jj) = BLK%DVLY(i,j,k)
          DVLZ_ALN(ii,jj) = BLK%DVLZ(i,j,k)
          !--------------------------------
          DWLX_ALN(ii,jj) = BLK%DWLX(i,j,k)
          DWLY_ALN(ii,jj) = BLK%DWLY(i,j,k)
          DWLZ_ALN(ii,jj) = BLK%DWLZ(i,j,k)
          !--------------------------------
          DTMX_ALN(ii,jj) = BLK%DTMX(i,j,k)
          DTMY_ALN(ii,jj) = BLK%DTMY(i,j,k)
          DTMZ_ALN(ii,jj) = BLK%DTMZ(i,j,k)
          !--------------------------------
          VSSV_ALN(ii,jj) = BLK%VSSV(i,j,k)
          VSST_ALN(ii,jj) = BLK%VSST(i,j,k)
          VSBV_ALN(ii,jj) = BLK%VSBV(i,j,k)
          VSBT_ALN(ii,jj) = BLK%VSBT(i,j,k)
          CKPV_ALN(ii,jj) = BLK%CKPV(i,j,k)
          CKPT_ALN(ii,jj) = BLK%CKPT(i,j,k)
          !--------------------------------
          DXSX_ALN(ii,jj) = BLK%DZTX(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DZTY(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DZTZ(i,j,k)
          !--------------------------------
          TMPR_ALN(ii,jj) = BLK%TMPR(i,j,k)  !! FURUSAWA
        enddo
      enddo

      do l=1,neqns-5
        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + (ii-1)/imax+1
            k = koft + jj
            YSPC_ALN(ii,jj,l) = BLK%YSPC(i,j,k,l)
            !------------------------------------
            DYSX_ALN(ii,jj,l) = BLK%DYSX(i,j,k,l)
            DYSY_ALN(ii,jj,l) = BLK%DYSY(i,j,k,l)
            DYSZ_ALN(ii,jj,l) = BLK%DYSZ(i,j,k,l)
          enddo
        enddo
      enddo


      call calc_viscous


      do l=1,neqns
        if (    (l==2).or.(l==3).or.(l==4).or.(l==5) &
            .or.(iact_tbl_sa>0.and.l==neqn_tbl_sa)   &
            .or.(l==neqn_passive_scalar_1)           &
            .or.(l==neqn_passive_scalar_2)           &
            .or.(l==neqn_passive_scalar_3)            ) then

          call calc_fdm_1st_align(FLUX_ALN(:,:,l),FLUX_ALN(:,:,l),0,icompact,isc_zta,iimax,jjmax)

        endif
      enddo


      do l=1,neqns
        if (    (l==2).or.(l==3).or.(l==4).or.(l==5) &
            .or.(iact_tbl_sa>0.and.l==neqn_tbl_sa)   &
            .or.(l==neqn_passive_scalar_1)           &
            .or.(l==neqn_passive_scalar_2)           &
            .or.(l==neqn_passive_scalar_3)            ) then

          do jj=1,kmax
            do ii=1,imax*jmax
              i = BLK%ista-1 + mod(ii-1,imax)+1
              j = BLK%jsta-1 + (ii-1)/imax+1
              k = BLK%ksta-1 + jj
              RHS(i,j,k,l) = RHS(i,j,k,l) - FLUX_ALN(ii,nhalo1+jj,l)
            enddo
          enddo

        endif
      enddo

      deallocate(DNST_ALN)
      deallocate(ULCT_ALN)
      deallocate(VLCT_ALN)
      deallocate(WLCT_ALN)
      if (neqns-5>0) then
        deallocate(YSPC_ALN)
      endif
      !-------------------
      deallocate(DULX_ALN)
      deallocate(DULY_ALN)
      deallocate(DULZ_ALN)
      !-------------------
      deallocate(DVLX_ALN)
      deallocate(DVLY_ALN)
      deallocate(DVLZ_ALN)
      !-------------------
      deallocate(DWLX_ALN)
      deallocate(DWLY_ALN)
      deallocate(DWLZ_ALN)
      !-------------------
      deallocate(DTMX_ALN)
      deallocate(DTMY_ALN)
      deallocate(DTMZ_ALN)
      !-------------------
      if (neqns-5>0) then
        deallocate(DYSX_ALN)
        deallocate(DYSY_ALN)
        deallocate(DYSZ_ALN)
      endif
      !-------------------
      deallocate(VSSV_ALN)
      deallocate(VSST_ALN)
      deallocate(VSBV_ALN)
      deallocate(VSBT_ALN)
      deallocate(CKPV_ALN)
      deallocate(CKPT_ALN)
      !-------------------
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      !-------------------
      deallocate(FLUX_ALN)
      !-------------------
      deallocate(TMPR_ALN)
    end subroutine calc_zta
!-------------------------------------------------------------------------------
    subroutine calc_viscous
      integer :: ii,jj

      real(8) :: rh,u,v,w
      real(8) :: ux,uy,uz,vx,vy,vz,wx,wy,wz,tmx,tmy,tmz
      real(8) :: xsix,xsiy,xsiz
      real(8) :: txx,tyy,tzz,txy,tyz,tzx,tt
      real(8) :: cmu,cbk,ckp
      real(8) :: utx,uty,utz
      real(8) :: rre

      real(8) :: tkx,tky,tkz
      real(8) :: cfftk

      real(8) :: cdf,cdft
      real(8) :: kex1,key1,kez1
      real(8) :: kex2,key2,kez2
      real(8) :: kex3,key3,kez3

      real(8) :: t, PDIFFI,schmt !! FURUSAWA


      do jj=1,jjmax
        do ii=1,iimax
          rh = DNST_ALN(ii,jj)
          u  = ULCT_ALN(ii,jj)
          v  = VLCT_ALN(ii,jj)
          w  = WLCT_ALN(ii,jj)
          !-------------------
          ux = DULX_ALN(ii,jj)
          uy = DULY_ALN(ii,jj)
          uz = DULZ_ALN(ii,jj)
          vx = DVLX_ALN(ii,jj)
          vy = DVLY_ALN(ii,jj)
          vz = DVLZ_ALN(ii,jj)
          wx = DWLX_ALN(ii,jj)
          wy = DWLY_ALN(ii,jj)
          wz = DWLZ_ALN(ii,jj)
          tmx= DTMX_ALN(ii,jj)
          tmy= DTMY_ALN(ii,jj)
          tmz= DTMZ_ALN(ii,jj)
          !-----------------------------------
          cmu= VSSV_ALN(ii,jj)+VSST_ALN(ii,jj)
          cbk= VSBV_ALN(ii,jj)+VSBT_ALN(ii,jj)
          ckp= CKPV_ALN(ii,jj)+CKPT_ALN(ii,jj)
          !-----------------------------------
          xsix = DXSX_ALN(ii,jj)
          xsiy = DXSY_ALN(ii,jj)
          xsiz = DXSZ_ALN(ii,jj)

          rre= 1.0d0/RYNLD

          ckp= ckp*CPREF

          tt = (cbk-2.0d0*cmu/3.0d0)*(ux+vy+wz)

          txx = 2.0d0*cmu*ux+tt
          tyy = 2.0d0*cmu*vy+tt
          tzz = 2.0d0*cmu*wz+tt

          txy = cmu*(uy+vx)
          tyz = cmu*(vz+wy)
          tzx = cmu*(wx+uz)

          utx= u*txx+v*txy+w*tzx+ckp*tmx
          uty= u*txy+v*tyy+w*tyz+ckp*tmy
          utz= u*tzx+v*tyz+w*tzz+ckp*tmz


          FLUX_ALN(ii,jj,2) = -rre*(xsix*txx+xsiy*txy+xsiz*tzx)
          FLUX_ALN(ii,jj,3) = -rre*(xsix*txy+xsiy*tyy+xsiz*tyz)
          FLUX_ALN(ii,jj,4) = -rre*(xsix*tzx+xsiy*tyz+xsiz*tzz)
          FLUX_ALN(ii,jj,5) = -rre*(xsix*utx+xsiy*uty+xsiz*utz)

          !Turbulence model
          !S-A
          if (iact_tbl_sa>0) then
            tkx = DYSX_ALN(ii,jj,neqn_tbl_sa-5)
            tky = DYSY_ALN(ii,jj,neqn_tbl_sa-5)
            tkz = DYSZ_ALN(ii,jj,neqn_tbl_sa-5)
            cfftk = (VSSV_ALN(ii,jj)*VSREF+rh*RHREF*YSPC_ALN(ii,jj,neqn_tbl_sa-5))/(SGMSA*VSREF)
            FLUX_ALN(ii,jj,neqn_tbl_sa) = -rre*cfftk*(xsix*tkx+xsiy*tky+xsiz*tkz)
          endif

          !Scalar diffusion

          t  =  TMPR_ALN(ii,jj)*TTREF !! FURUSAWA
!!F          PDIFFI = 10.0D0**(2.67D0)*DEXP(-15.9D3/(8.3144D0*t))*1.0D-9
          PDIFFI = 2.24D-13*(t)**0.763/(rh*RHREF)
          schmt = VSSV_ALN(ii,jj)*VSREF/(rh*RHREF *PDIFFI)  !! FURUSAWA
          
          cdf = VSSV_ALN(ii,jj)*VSREF/(rh*RHREF* schmt)/DFREF
          cdft= VSST_ALN(ii,jj)*VSREF/(rh*RHREF*tschmt)/DFREF
          if (neqn_passive_scalar_1>0) then
            kex1 = (cdf+cdft)*rh*DYSX_ALN(ii,jj,neqn_passive_scalar_1-5)
            key1 = (cdf+cdft)*rh*DYSY_ALN(ii,jj,neqn_passive_scalar_1-5)
            kez1 = (cdf+cdft)*rh*DYSZ_ALN(ii,jj,neqn_passive_scalar_1-5)
            FLUX_ALN(ii,jj,neqn_passive_scalar_1) = -1.0d0*(xsix*kex1+xsiy*key1+xsiz*kez1)
          endif
          if (neqn_passive_scalar_2>0) then
            kex2 = (cdf+cdft)*rh*DYSX_ALN(ii,jj,neqn_passive_scalar_2-5)
            key2 = (cdf+cdft)*rh*DYSY_ALN(ii,jj,neqn_passive_scalar_2-5)
            kez2 = (cdf+cdft)*rh*DYSZ_ALN(ii,jj,neqn_passive_scalar_2-5)
            FLUX_ALN(ii,jj,neqn_passive_scalar_2) = -1.0d0*(xsix*kex2+xsiy*key2+xsiz*kez2)
          endif
          if (neqn_passive_scalar_3>0) then
            kex3 = (cdf+cdft)*rh*DYSX_ALN(ii,jj,neqn_passive_scalar_3-5)
            key3 = (cdf+cdft)*rh*DYSY_ALN(ii,jj,neqn_passive_scalar_3-5)
            kez3 = (cdf+cdft)*rh*DYSZ_ALN(ii,jj,neqn_passive_scalar_3-5)
            FLUX_ALN(ii,jj,neqn_passive_scalar_3) = -1.0d0*(xsix*kex3+xsiy*key3+xsiz*kez3)
          endif

        enddo
      enddo

    end subroutine calc_viscous
!-------------------------------------------------------------------------------
  end subroutine flux_viscous

end module mdl_flux_viscous
