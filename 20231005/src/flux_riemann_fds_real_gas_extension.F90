module mdl_flux_riemann_fds_real_gas_extension
  use mdl_mpisub_sbsp
  use mdl_decompo
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_muscl
  use mdl_weno
  use mdl_halo
  use mdl_table
  implicit none

!-------------------------------------------------------------------------------
!Flux difference splitting scheme with MUSCL or 5th-order WENO.
!Real gas extended version.
!-------------------------------------------------------------------------------

  integer, private, parameter :: isc_reconstruction = 0 !<1>:weno5 <else>:MUSCL

  integer, private, parameter :: iorder_muscl = 3

  integer, private, parameter :: iact_table_thermo = 0

contains

  subroutine flux_riemann_fds_real_gas_extension(BLK,RHS)
    type(block), intent(inout) :: BLK
    real(8), intent(inout), dimension(:,:,:,:) :: RHS

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: iimax,jjmax

    !Cell center
    real(8), allocatable, dimension(:,:,:) :: QQ_ALN !T,u,v,w,p,ys1,...

    !Cell face
    real(8), allocatable, dimension(:,:,:) :: QL_ALN
    real(8), allocatable, dimension(:,:,:) :: QR_ALN
    !-------------------------------------------------
    real(8), allocatable, dimension(:,:)   :: DXSX_ALN
    real(8), allocatable, dimension(:,:)   :: DXSY_ALN
    real(8), allocatable, dimension(:,:)   :: DXSZ_ALN
    !-------------------------------------------------
    real(8), allocatable, dimension(:,:)   :: XCRD_ALN
    real(8), allocatable, dimension(:,:)   :: YCRD_ALN
    real(8), allocatable, dimension(:,:)   :: ZCRD_ALN
    !-------------------------------------------------
    real(8), allocatable, dimension(:,:,:) :: FLUX_ALN
    real(8), allocatable, dimension(:,:,:) :: DLUX_ALN

    real(8), allocatable, dimension(:,:)   :: DNST_L_ALN
    real(8), allocatable, dimension(:,:)   :: DNST_R_ALN
    !---------------------------------------------------
    real(8), allocatable, dimension(:,:)   :: HTPS_L_ALN
    real(8), allocatable, dimension(:,:)   :: HTPS_R_ALN
    !---------------------------------------------------
    real(8), allocatable, dimension(:,:)   :: SHCP_L_ALN
    real(8), allocatable, dimension(:,:)   :: SHCP_R_ALN
    !---------------------------------------------------
    real(8), allocatable, dimension(:,:)   :: DHPR_L_ALN
    real(8), allocatable, dimension(:,:)   :: DHPR_R_ALN
    !---------------------------------------------------
    real(8), allocatable, dimension(:,:)   :: DRTM_L_ALN
    real(8), allocatable, dimension(:,:)   :: DRTM_R_ALN
    !---------------------------------------------------
    real(8), allocatable, dimension(:,:)   :: DRPR_L_ALN
    real(8), allocatable, dimension(:,:)   :: DRPR_R_ALN

    !Cell center
    real(8), allocatable, dimension(:,:)   :: DNST_ALN
    real(8), allocatable, dimension(:,:)   :: HTPS_ALN
    real(8), allocatable, dimension(:,:)   :: SHCP_ALN
    real(8), allocatable, dimension(:,:)   :: DHPR_ALN
    real(8), allocatable, dimension(:,:)   :: DRTM_ALN
    real(8), allocatable, dimension(:,:)   :: DRPR_ALN
    real(8), allocatable, dimension(:,:)   :: VPRE_ALN
    !-------------------------------------------------
    real(8), allocatable, dimension(:,:,:) :: DFLX_ALN

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call calc_xsi
    call calc_eta
    if (.not.(i2dimension>0)) then
      call calc_zta
    endif

  contains

    subroutine calc_xsi
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft
      integer :: ivalid

      nhalo1 = 0
      nhalo2 = 0

      if (IFLAG_PRIM(1)==1) nhalo1 = halo_width
      if (IFLAG_PRIM(2)==1) nhalo2 = halo_width

      iimax = jmax*kmax
      jjmax = imax+nhalo1+nhalo2

      ioft = BLK%ista-1
      joft = BLK%jsta-1
      koft = BLK%ksta-1

      if (IFLAG_PRIM(1)==1) ioft = 0

      allocate(QQ_ALN(iimax,jjmax,neqns))

      !Cell face----------------------------------
      !reconstructed values
      allocate(QL_ALN(iimax,jjmax+1,neqns))
      allocate(QR_ALN(iimax,jjmax+1,neqns))
      !xsix,xsiy,xsiz
      allocate(DXSX_ALN(iimax,jjmax+1))
      allocate(DXSY_ALN(iimax,jjmax+1))
      allocate(DXSZ_ALN(iimax,jjmax+1))
      !x,y,z coordinate
      allocate(XCRD_ALN(iimax,jjmax+1))
      allocate(YCRD_ALN(iimax,jjmax+1))
      allocate(ZCRD_ALN(iimax,jjmax+1))
      !advection flux
      allocate(FLUX_ALN(iimax,jjmax+1,neqns))
      !dissipation flux
      allocate(DLUX_ALN(iimax,jjmax+1,neqns))
      !density
      allocate(DNST_L_ALN(iimax,jjmax+1))
      allocate(DNST_R_ALN(iimax,jjmax+1))
      !static enthalpy
      allocate(HTPS_L_ALN(iimax,jjmax+1))
      allocate(HTPS_R_ALN(iimax,jjmax+1))

      if (iact_table_thermo==1) then
        !Cp
        allocate(SHCP_L_ALN(iimax,jjmax+1))
        allocate(SHCP_R_ALN(iimax,jjmax+1))
        !h,p
        allocate(DHPR_L_ALN(iimax,jjmax+1))
        allocate(DHPR_R_ALN(iimax,jjmax+1))
        !rh,T
        allocate(DRTM_L_ALN(iimax,jjmax+1))
        allocate(DRTM_R_ALN(iimax,jjmax+1))
        !rh,p
        allocate(DRPR_L_ALN(iimax,jjmax+1))
        allocate(DRPR_R_ALN(iimax,jjmax+1))
      endif

      !Cell center--------------------------------
      !density
      allocate(DNST_ALN(iimax,jjmax))
      !static enthalpy
      allocate(HTPS_ALN(iimax,jjmax))
      !Cp
      allocate(SHCP_ALN(iimax,jjmax))
      !h,p
      allocate(DHPR_ALN(iimax,jjmax))
      !rh,T
      allocate(DRTM_ALN(iimax,jjmax))
      !rh,p
      allocate(DRPR_ALN(iimax,jjmax))
      !Velocity scale for preconditioning
      allocate(VPRE_ALN(iimax,jjmax))
      !rhs
      allocate(DFLX_ALN(iimax,jjmax,neqns))

      do jj=1,jjmax+1
        do ii=1,iimax
          i = ioft + jj
          j = joft + mod(ii-1,jmax)+1
          k = koft + (ii-1)/jmax+1
          DXSX_ALN(ii,jj) = BLK%DXSX_FACE(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DXSY_FACE(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DXSZ_FACE(i,j,k)
          !--------------------------------
          XCRD_ALN(ii,jj) = BLK%XCRD_FACE_ETAZTA(i,j,k)
          YCRD_ALN(ii,jj) = BLK%YCRD_FACE_ETAZTA(i,j,k)
          ZCRD_ALN(ii,jj) = BLK%ZCRD_FACE_ETAZTA(i,j,k)
        enddo
      enddo

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + jj
          j = joft + mod(ii-1,jmax)+1
          k = koft + (ii-1)/jmax+1
          QQ_ALN(ii,jj,1) = BLK%TMPR(i,j,k)
          QQ_ALN(ii,jj,2) = BLK%ULCT(i,j,k)
          QQ_ALN(ii,jj,3) = BLK%VLCT(i,j,k)
          QQ_ALN(ii,jj,4) = BLK%WLCT(i,j,k)
          QQ_ALN(ii,jj,5) = BLK%PRSS(i,j,k)
          !
          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          HTPS_ALN(ii,jj) = BLK%HTPS(i,j,k)
          SHCP_ALN(ii,jj) = BLK%SHCP(i,j,k)
          DHPR_ALN(ii,jj) = BLK%DHPR(i,j,k)
          DRTM_ALN(ii,jj) = BLK%DRTM(i,j,k)
          DRPR_ALN(ii,jj) = BLK%DRPR(i,j,k)
          VPRE_ALN(ii,jj) = BLK%VPRE(i,j,k)
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

      if (isc_reconstruction==1) then
        do l=1,neqns
          call weno5_align(QQ_ALN(:,:,l),QL_ALN(:,:,l),QR_ALN(:,:,l),iimax,jjmax)
        enddo
      else
        do l=1,neqns
          call muscl_align(QQ_ALN(:,:,l),QL_ALN(:,:,l),QR_ALN(:,:,l),iorder_muscl,iimax,jjmax)
        enddo
      endif

      !density
      call look_up_table_calc_func_2d(TBL_TP_to_DNST,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                      DNST_L_ALN(:,:),RHREF,1,iimax,1,jjmax+1,ivalid)
      call look_up_table_calc_func_2d(TBL_TP_to_DNST,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                      DNST_R_ALN(:,:),RHREF,1,iimax,1,jjmax+1,ivalid)

      !static enthalpy
      call look_up_table_calc_func_2d(TBL_TP_to_HTPS,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                      HTPS_L_ALN(:,:),TTREF,1,iimax,1,jjmax+1,ivalid)
      call look_up_table_calc_func_2d(TBL_TP_to_HTPS,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                      HTPS_R_ALN(:,:),TTREF,1,iimax,1,jjmax+1,ivalid)

      if (iact_table_thermo==1) then
        !Cp
        call look_up_table_calc_func_2d(TBL_TP_to_SHCP,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        SHCP_L_ALN(:,:),CPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_SHCP,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        SHCP_R_ALN(:,:),CPREF,1,iimax,1,jjmax+1,ivalid)

        !h,p
        call look_up_table_calc_func_2d(TBL_TP_to_DHPR,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DHPR_L_ALN(:,:),TTREF/PPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DHPR,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DHPR_R_ALN(:,:),TTREF/PPREF,1,iimax,1,jjmax+1,ivalid)

        !rh,T
        call look_up_table_calc_func_2d(TBL_TP_to_DRTM,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DRTM_L_ALN(:,:),RHREF/TTREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DRTM,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DRTM_R_ALN(:,:),RHREF/TTREF,1,iimax,1,jjmax+1,ivalid)

        !rh,p
        call look_up_table_calc_func_2d(TBL_TP_to_DRPR,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DRPR_L_ALN(:,:),RHREF/PPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DRPR,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DRPR_R_ALN(:,:),RHREF/PPREF,1,iimax,1,jjmax+1,ivalid)
      endif

      call calc_fds(1.0d0,0.0d0,0.0d0)

      do l=1,neqns
        do jj=1,imax
          do ii=1,jmax*kmax
            i = BLK%ista-1 + jj
            j = BLK%jsta-1 + mod(ii-1,jmax)+1
            k = BLK%ksta-1 + (ii-1)/jmax+1
            RHS(i,j,k,l) = RHS(i,j,k,l) + DFLX_ALN(ii,nhalo1+jj,l)
          enddo
        enddo
      enddo

      deallocate(QQ_ALN)
      deallocate(QL_ALN)
      deallocate(QR_ALN)
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      deallocate(XCRD_ALN)
      deallocate(YCRD_ALN)
      deallocate(ZCRD_ALN)
      deallocate(FLUX_ALN)
      deallocate(DLUX_ALN)

      deallocate(DNST_L_ALN)
      deallocate(DNST_R_ALN)
      deallocate(HTPS_L_ALN)
      deallocate(HTPS_R_ALN)

      if (iact_table_thermo==1) then
        deallocate(SHCP_L_ALN)
        deallocate(SHCP_R_ALN)
        deallocate(DHPR_L_ALN)
        deallocate(DHPR_R_ALN)
        deallocate(DRTM_L_ALN)
        deallocate(DRTM_R_ALN)
        deallocate(DRPR_L_ALN)
        deallocate(DRPR_R_ALN)
      endif

      deallocate(DNST_ALN)
      deallocate(HTPS_ALN)
      deallocate(SHCP_ALN)
      deallocate(DHPR_ALN)
      deallocate(DRTM_ALN)
      deallocate(DRPR_ALN)
      deallocate(VPRE_ALN)
      deallocate(DFLX_ALN)
    end subroutine calc_xsi
!-------------------------------------------------------------------------------
    subroutine calc_eta
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft
      integer :: ivalid

      nhalo1 = 0
      nhalo2 = 0

      if (IFLAG_PRIM(3)==1) nhalo1 = halo_width
      if (IFLAG_PRIM(4)==1) nhalo2 = halo_width

      iimax = imax*kmax
      jjmax = jmax+nhalo1+nhalo2

      ioft = BLK%ista-1
      joft = BLK%jsta-1
      koft = BLK%ksta-1

      if (IFLAG_PRIM(3)==1) joft = 0

      allocate(QQ_ALN(iimax,jjmax,neqns))

      !Cell face----------------------------------
      !reconstructed values
      allocate(QL_ALN(iimax,jjmax+1,neqns))
      allocate(QR_ALN(iimax,jjmax+1,neqns))
      !xsix,xsiy,xsiz
      allocate(DXSX_ALN(iimax,jjmax+1))
      allocate(DXSY_ALN(iimax,jjmax+1))
      allocate(DXSZ_ALN(iimax,jjmax+1))
      !x,y,z coordinate
      allocate(XCRD_ALN(iimax,jjmax+1))
      allocate(YCRD_ALN(iimax,jjmax+1))
      allocate(ZCRD_ALN(iimax,jjmax+1))
      !advection flux
      allocate(FLUX_ALN(iimax,jjmax+1,neqns))
      !dissipation flux
      allocate(DLUX_ALN(iimax,jjmax+1,neqns))
      !density
      allocate(DNST_L_ALN(iimax,jjmax+1))
      allocate(DNST_R_ALN(iimax,jjmax+1))
      !static enthalpy
      allocate(HTPS_L_ALN(iimax,jjmax+1))
      allocate(HTPS_R_ALN(iimax,jjmax+1))

      if (iact_table_thermo==1) then
        !Cp
        allocate(SHCP_L_ALN(iimax,jjmax+1))
        allocate(SHCP_R_ALN(iimax,jjmax+1))
        !h,p
        allocate(DHPR_L_ALN(iimax,jjmax+1))
        allocate(DHPR_R_ALN(iimax,jjmax+1))
        !rh,T
        allocate(DRTM_L_ALN(iimax,jjmax+1))
        allocate(DRTM_R_ALN(iimax,jjmax+1))
        !rh,p
        allocate(DRPR_L_ALN(iimax,jjmax+1))
        allocate(DRPR_R_ALN(iimax,jjmax+1))
      endif

      !Cell center--------------------------------
      !density
      allocate(DNST_ALN(iimax,jjmax))
      !static enthalpy
      allocate(HTPS_ALN(iimax,jjmax))
      !Cp
      allocate(SHCP_ALN(iimax,jjmax))
      !h,p
      allocate(DHPR_ALN(iimax,jjmax))
      !rh,T
      allocate(DRTM_ALN(iimax,jjmax))
      !rh,p
      allocate(DRPR_ALN(iimax,jjmax))
      !Velocity scale for preconditioning
      allocate(VPRE_ALN(iimax,jjmax))
      !rhs
      allocate(DFLX_ALN(iimax,jjmax,neqns))

      do jj=1,jjmax+1
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + jj
          k = koft + (ii-1)/imax+1
          DXSX_ALN(ii,jj) = BLK%DETX_FACE(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DETY_FACE(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DETZ_FACE(i,j,k)
          !--------------------------------
          XCRD_ALN(ii,jj) = BLK%XCRD_FACE_ZTAXSI(i,j,k)
          YCRD_ALN(ii,jj) = BLK%YCRD_FACE_ZTAXSI(i,j,k)
          ZCRD_ALN(ii,jj) = BLK%ZCRD_FACE_ZTAXSI(i,j,k)
        enddo
      enddo

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + jj
          k = koft + (ii-1)/imax+1
          QQ_ALN(ii,jj,1) = BLK%TMPR(i,j,k)
          QQ_ALN(ii,jj,2) = BLK%ULCT(i,j,k)
          QQ_ALN(ii,jj,3) = BLK%VLCT(i,j,k)
          QQ_ALN(ii,jj,4) = BLK%WLCT(i,j,k)
          QQ_ALN(ii,jj,5) = BLK%PRSS(i,j,k)
          !
          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          HTPS_ALN(ii,jj) = BLK%HTPS(i,j,k)
          SHCP_ALN(ii,jj) = BLK%SHCP(i,j,k)
          DHPR_ALN(ii,jj) = BLK%DHPR(i,j,k)
          DRTM_ALN(ii,jj) = BLK%DRTM(i,j,k)
          DRPR_ALN(ii,jj) = BLK%DRPR(i,j,k)
          VPRE_ALN(ii,jj) = BLK%VPRE(i,j,k)
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

      if (isc_reconstruction==1) then
        do l=1,neqns
          call weno5_align(QQ_ALN(:,:,l),QL_ALN(:,:,l),QR_ALN(:,:,l),iimax,jjmax)
        enddo
      else
        do l=1,neqns
          call muscl_align(QQ_ALN(:,:,l),QL_ALN(:,:,l),QR_ALN(:,:,l),iorder_muscl,iimax,jjmax)
        enddo
      endif

      !density
      call look_up_table_calc_func_2d(TBL_TP_to_DNST,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                      DNST_L_ALN(:,:),RHREF,1,iimax,1,jjmax+1,ivalid)
      call look_up_table_calc_func_2d(TBL_TP_to_DNST,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                      DNST_R_ALN(:,:),RHREF,1,iimax,1,jjmax+1,ivalid)

      !static enthalpy
      call look_up_table_calc_func_2d(TBL_TP_to_HTPS,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                      HTPS_L_ALN(:,:),TTREF,1,iimax,1,jjmax+1,ivalid)
      call look_up_table_calc_func_2d(TBL_TP_to_HTPS,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                      HTPS_R_ALN(:,:),TTREF,1,iimax,1,jjmax+1,ivalid)

      if (iact_table_thermo==1) then
        !Cp
        call look_up_table_calc_func_2d(TBL_TP_to_SHCP,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        SHCP_L_ALN(:,:),CPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_SHCP,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        SHCP_R_ALN(:,:),CPREF,1,iimax,1,jjmax+1,ivalid)

        !h,p
        call look_up_table_calc_func_2d(TBL_TP_to_DHPR,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DHPR_L_ALN(:,:),TTREF/PPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DHPR,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DHPR_R_ALN(:,:),TTREF/PPREF,1,iimax,1,jjmax+1,ivalid)

        !rh,T
        call look_up_table_calc_func_2d(TBL_TP_to_DRTM,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DRTM_L_ALN(:,:),RHREF/TTREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DRTM,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DRTM_R_ALN(:,:),RHREF/TTREF,1,iimax,1,jjmax+1,ivalid)

        !rh,p
        call look_up_table_calc_func_2d(TBL_TP_to_DRPR,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DRPR_L_ALN(:,:),RHREF/PPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DRPR,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DRPR_R_ALN(:,:),RHREF/PPREF,1,iimax,1,jjmax+1,ivalid)
      endif

      call calc_fds(0.0d0,1.0d0,0.0d0)

      do l=1,neqns
        do jj=1,jmax
          do ii=1,imax*kmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + jj
            k = BLK%ksta-1 + (ii-1)/imax+1
            RHS(i,j,k,l) = RHS(i,j,k,l) + DFLX_ALN(ii,nhalo1+jj,l)
          enddo
        enddo
      enddo

      deallocate(QQ_ALN)
      deallocate(QL_ALN)
      deallocate(QR_ALN)
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      deallocate(XCRD_ALN)
      deallocate(YCRD_ALN)
      deallocate(ZCRD_ALN)
      deallocate(FLUX_ALN)
      deallocate(DLUX_ALN)

      deallocate(DNST_L_ALN)
      deallocate(DNST_R_ALN)
      deallocate(HTPS_L_ALN)
      deallocate(HTPS_R_ALN)

      if (iact_table_thermo==1) then
        deallocate(SHCP_L_ALN)
        deallocate(SHCP_R_ALN)
        deallocate(DHPR_L_ALN)
        deallocate(DHPR_R_ALN)
        deallocate(DRTM_L_ALN)
        deallocate(DRTM_R_ALN)
        deallocate(DRPR_L_ALN)
        deallocate(DRPR_R_ALN)
      endif

      deallocate(DNST_ALN)
      deallocate(HTPS_ALN)
      deallocate(SHCP_ALN)
      deallocate(DHPR_ALN)
      deallocate(DRTM_ALN)
      deallocate(DRPR_ALN)
      deallocate(VPRE_ALN)
      deallocate(VPRE_ALN)
      deallocate(DFLX_ALN)
    end subroutine calc_eta
!-------------------------------------------------------------------------------
    subroutine calc_zta
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft
      integer :: ivalid

      nhalo1 = 0
      nhalo2 = 0

      if (IFLAG_PRIM(5)==1) nhalo1 = halo_width
      if (IFLAG_PRIM(6)==1) nhalo2 = halo_width

      iimax = imax*jmax
      jjmax = kmax+nhalo1+nhalo2

      ioft = BLK%ista-1
      joft = BLK%jsta-1
      koft = BLK%ksta-1

      if (IFLAG_PRIM(5)==1) koft = 0

      allocate(QQ_ALN(iimax,jjmax,neqns))

      !Cell face----------------------------------
      !reconstructed values
      allocate(QL_ALN(iimax,jjmax+1,neqns))
      allocate(QR_ALN(iimax,jjmax+1,neqns))
      !xsix,xsiy,xsiz
      allocate(DXSX_ALN(iimax,jjmax+1))
      allocate(DXSY_ALN(iimax,jjmax+1))
      allocate(DXSZ_ALN(iimax,jjmax+1))
      !x,y,z coordinate
      allocate(XCRD_ALN(iimax,jjmax+1))
      allocate(YCRD_ALN(iimax,jjmax+1))
      allocate(ZCRD_ALN(iimax,jjmax+1))
      !advection flux
      allocate(FLUX_ALN(iimax,jjmax+1,neqns))
      !dissipation flux
      allocate(DLUX_ALN(iimax,jjmax+1,neqns))
      !density
      allocate(DNST_L_ALN(iimax,jjmax+1))
      allocate(DNST_R_ALN(iimax,jjmax+1))
      !static enthalpy
      allocate(HTPS_L_ALN(iimax,jjmax+1))
      allocate(HTPS_R_ALN(iimax,jjmax+1))

      if (iact_table_thermo==1) then
        !Cp
        allocate(SHCP_L_ALN(iimax,jjmax+1))
        allocate(SHCP_R_ALN(iimax,jjmax+1))
        !h,p
        allocate(DHPR_L_ALN(iimax,jjmax+1))
        allocate(DHPR_R_ALN(iimax,jjmax+1))
        !rh,T
        allocate(DRTM_L_ALN(iimax,jjmax+1))
        allocate(DRTM_R_ALN(iimax,jjmax+1))
        !rh,p
        allocate(DRPR_L_ALN(iimax,jjmax+1))
        allocate(DRPR_R_ALN(iimax,jjmax+1))
      endif

      !Cell center--------------------------------
      !density
      allocate(DNST_ALN(iimax,jjmax))
      !static enthalpy
      allocate(HTPS_ALN(iimax,jjmax))
      !Cp
      allocate(SHCP_ALN(iimax,jjmax))
      !h,p
      allocate(DHPR_ALN(iimax,jjmax))
      !rh,T
      allocate(DRTM_ALN(iimax,jjmax))
      !rh,p
      allocate(DRPR_ALN(iimax,jjmax))
      !Velocity scale for preconditioning
      allocate(VPRE_ALN(iimax,jjmax))
      !rhs
      allocate(DFLX_ALN(iimax,jjmax,neqns))

      do jj=1,jjmax+1
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + (ii-1)/imax+1
          k = koft + jj
          DXSX_ALN(ii,jj) = BLK%DZTX_FACE(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DZTY_FACE(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DZTZ_FACE(i,j,k)
          !--------------------------------
          XCRD_ALN(ii,jj) = BLK%XCRD_FACE_XSIETA(i,j,k)
          YCRD_ALN(ii,jj) = BLK%YCRD_FACE_XSIETA(i,j,k)
          ZCRD_ALN(ii,jj) = BLK%ZCRD_FACE_XSIETA(i,j,k)
        enddo
      enddo

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + (ii-1)/imax+1
          k = koft + jj
          QQ_ALN(ii,jj,1) = BLK%TMPR(i,j,k)
          QQ_ALN(ii,jj,2) = BLK%ULCT(i,j,k)
          QQ_ALN(ii,jj,3) = BLK%VLCT(i,j,k)
          QQ_ALN(ii,jj,4) = BLK%WLCT(i,j,k)
          QQ_ALN(ii,jj,5) = BLK%PRSS(i,j,k)
          !
          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          HTPS_ALN(ii,jj) = BLK%HTPS(i,j,k)
          SHCP_ALN(ii,jj) = BLK%SHCP(i,j,k)
          DHPR_ALN(ii,jj) = BLK%DHPR(i,j,k)
          DRTM_ALN(ii,jj) = BLK%DRTM(i,j,k)
          DRPR_ALN(ii,jj) = BLK%DRPR(i,j,k)
          VPRE_ALN(ii,jj) = BLK%VPRE(i,j,k)
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

      if (isc_reconstruction==1) then
        do l=1,neqns
          call weno5_align(QQ_ALN(:,:,l),QL_ALN(:,:,l),QR_ALN(:,:,l),iimax,jjmax)
        enddo
      else
        do l=1,neqns
          call muscl_align(QQ_ALN(:,:,l),QL_ALN(:,:,l),QR_ALN(:,:,l),iorder_muscl,iimax,jjmax)
        enddo
      endif

      !density
      call look_up_table_calc_func_2d(TBL_TP_to_DNST,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                      DNST_L_ALN(:,:),RHREF,1,iimax,1,jjmax+1,ivalid)
      call look_up_table_calc_func_2d(TBL_TP_to_DNST,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                      DNST_R_ALN(:,:),RHREF,1,iimax,1,jjmax+1,ivalid)

      !static enthalpy
      call look_up_table_calc_func_2d(TBL_TP_to_HTPS,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                      HTPS_L_ALN(:,:),TTREF,1,iimax,1,jjmax+1,ivalid)
      call look_up_table_calc_func_2d(TBL_TP_to_HTPS,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                      HTPS_R_ALN(:,:),TTREF,1,iimax,1,jjmax+1,ivalid)

      if (iact_table_thermo==1) then
        !Cp
        call look_up_table_calc_func_2d(TBL_TP_to_SHCP,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        SHCP_L_ALN(:,:),CPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_SHCP,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        SHCP_R_ALN(:,:),CPREF,1,iimax,1,jjmax+1,ivalid)

        !h,p
        call look_up_table_calc_func_2d(TBL_TP_to_DHPR,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DHPR_L_ALN(:,:),TTREF/PPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DHPR,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DHPR_R_ALN(:,:),TTREF/PPREF,1,iimax,1,jjmax+1,ivalid)

        !rh,T
        call look_up_table_calc_func_2d(TBL_TP_to_DRTM,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DRTM_L_ALN(:,:),RHREF/TTREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DRTM,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DRTM_R_ALN(:,:),RHREF/TTREF,1,iimax,1,jjmax+1,ivalid)

        !rh,p
        call look_up_table_calc_func_2d(TBL_TP_to_DRPR,QL_ALN(:,:,1),TTREF,QL_ALN(:,:,5),PPREF, &
                                        DRPR_L_ALN(:,:),RHREF/PPREF,1,iimax,1,jjmax+1,ivalid)
        call look_up_table_calc_func_2d(TBL_TP_to_DRPR,QR_ALN(:,:,1),TTREF,QR_ALN(:,:,5),PPREF, &
                                        DRPR_R_ALN(:,:),RHREF/PPREF,1,iimax,1,jjmax+1,ivalid)
      endif

      call calc_fds(0.0d0,0.0d0,1.0d0)

      do l=1,neqns
        do jj=1,kmax
          do ii=1,imax*jmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + (ii-1)/imax+1
            k = BLK%ksta-1 + jj
            RHS(i,j,k,l) = RHS(i,j,k,l) + DFLX_ALN(ii,nhalo1+jj,l)
          enddo
        enddo
      enddo

      deallocate(QQ_ALN)
      deallocate(QL_ALN)
      deallocate(QR_ALN)
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      deallocate(XCRD_ALN)
      deallocate(YCRD_ALN)
      deallocate(ZCRD_ALN)
      deallocate(FLUX_ALN)
      deallocate(DLUX_ALN)

      deallocate(DNST_L_ALN)
      deallocate(DNST_R_ALN)
      deallocate(HTPS_L_ALN)
      deallocate(HTPS_R_ALN)

      if (iact_table_thermo==1) then
        deallocate(SHCP_L_ALN)
        deallocate(SHCP_R_ALN)
        deallocate(DHPR_L_ALN)
        deallocate(DHPR_R_ALN)
        deallocate(DRTM_L_ALN)
        deallocate(DRTM_R_ALN)
        deallocate(DRPR_L_ALN)
        deallocate(DRPR_R_ALN)
      endif

      deallocate(DNST_ALN)
      deallocate(HTPS_ALN)
      deallocate(SHCP_ALN)
      deallocate(DHPR_ALN)
      deallocate(DRTM_ALN)
      deallocate(DRPR_ALN)
      deallocate(VPRE_ALN)
      deallocate(VPRE_ALN)
      deallocate(DFLX_ALN)
    end subroutine calc_zta
!-------------------------------------------------------------------------------
    subroutine calc_fds(cd1,cd2,cd3)
      real(8), intent(in) :: cd1,cd2,cd3

      integer :: ii,jj,jl,jr,l
      real(8) :: sxm,sym,szm,sxyzm2,sxyzm,smm
      real(8) :: rhl,ul,vl,wl,hsl,htl,qql,rhtl,rhpl,hstl,hspl,tl,pl,ucl,c2l,cl,vpl,vp2l,ggdl
      real(8) :: rhr,ur,vr,wr,hsr,htr,qqr,rhtr,rhpr,hstr,hspr,tr,pr,ucr,c2r,cr,vpr,vp2r,ggdr
      real(8) :: ys01l,ys02l,ys03l,ys04l,ys05l,ys06l,ys07l
      real(8) :: ys01r,ys02r,ys03r,ys04r,ys05r,ys06r,ys07r
      real(8) :: f01l,f02l,f03l,f04l,f05l,f06l,f07l,f08l,f09l,f10l,f11l,f12l
      real(8) :: f01r,f02r,f03r,f04r,f05r,f06r,f07r,f08r,f09r,f10r,f11r,f12r
      real(8) :: g01l,g02l,g03l,g04l,g05l,g06l,g07l,g08l,g09l,g10l,g11l,g12l
      real(8) :: g01r,g02r,g03r,g04r,g05r,g06r,g07r,g08r,g09r,g10r,g11r,g12r
      real(8) :: q01l,q02l,q03l,q04l,q05l,q06l,q07l,q08l,q09l,q10l,q11l,q12l
      real(8) :: q01r,q02r,q03r,q04r,q05r,q06r,q07r,q08r,q09r,q10r,q11r,q12r
      real(8) :: s01l,s02l,s03l,s04l,s05l,s06l,s07l,s08l,s09l,s10l,s11l,s12l
      real(8) :: s01r,s02r,s03r,s04r,s05r,s06r,s07r,s08r,s09r,s10r,s11r,s12r
      real(8) :: eigal,eigbl,eigcl
      real(8) :: eigar,eigbr,eigcr
      real(8) :: eig1l,eig2l,eig3l,eig4l,eig5l
      real(8) :: eig1r,eig2r,eig3r,eig4r,eig5r
      real(8) :: eig1max,eig2max,eig3max,eig4max,eig5max
      real(8) :: llf01l,llf02l,llf03l,llf04l,llf05l,llf06l,llf07l,llf08l,llf09l,llf10l,llf11l,llf12l
      real(8) :: llf01r,llf02r,llf03r,llf04r,llf05r,llf06r,llf07r,llf08r,llf09r,llf10r,llf11r,llf12r
      real(8) :: sig1,sig2,sig3,sig4,sig5
      real(8) :: sig1l,sig2l,sig3l,sig4l,sig5l
      real(8) :: sig1r,sig2r,sig3r,sig4r,sig5r
      real(8) :: xai01l,xai02l,xai03l,xai04l,xai05l,xai06l,xai07l,xai08l,xai09l,xai10l,xai11l,xai12l
      real(8) :: xai01r,xai02r,xai03r,xai04r,xai05r,xai06r,xai07r,xai08r,xai09r,xai10r,xai11r,xai12r

      real(8) :: tht0,tht2,tht3,thtb,thtb0,thtb2,thtb3

      real(8) :: m11,m12,m13,m14,m15
      real(8) :: m21,m22,m23,m24,m25
      real(8) :: m31,m32,m33,m34,m35
      real(8) :: m41,m42,m43,m44,m45
      real(8) :: m51,m52,m53,m54,m55

      real(8) :: alph,rrh,hsqq,dp,rdp

      real(8) :: pi
      real(8) :: nx,ny,nz,ox,oy,oz
      real(8) :: am,bm,gm,acm,rxm,rym,rzm


      pi = 4.0d0*datan(1.0d0)

      nx = 0.0d0
      ny = 0.0d0
      nz = 0.0d0
      if (frame_axis==1) nx = 1.0d0
      if (frame_axis==2) ny = 1.0d0
      if (frame_axis==3) nz = 1.0d0

      ox = 2.0d0*pi*frame_rpm/60.0d0*nx*DTREF
      oy = 2.0d0*pi*frame_rpm/60.0d0*ny*DTREF
      oz = 2.0d0*pi*frame_rpm/60.0d0*nz*DTREF

      do jj=1,jjmax+1
        jl = max0(jj-1,    1)
        jr = min0(jj  ,jjmax)
        do ii=1,iimax
          sxm = DXSX_ALN(ii,jj)
          sym = DXSY_ALN(ii,jj)
          szm = DXSZ_ALN(ii,jj)
          sxyzm2 = sxm*sxm+sym*sym+szm*szm
          sxyzm  = dsqrt(sxyzm2)
          smm = cd1*sxm+cd2*sym+cd3*szm

          rxm = XCRD_ALN(ii,jj)-frame_origin(1)
          rym = YCRD_ALN(ii,jj)-frame_origin(2)
          rzm = ZCRD_ALN(ii,jj)-frame_origin(3)

          am = oy*rzm-oz*rym
          bm = oz*rxm-ox*rzm
          gm = ox*rym-oy*rxm
          acm= am*sxm+bm*sym+gm*szm

          !left side
          tl = QL_ALN(ii,jj,1)
          ul = QL_ALN(ii,jj,2)
          vl = QL_ALN(ii,jj,3)
          wl = QL_ALN(ii,jj,4)
          pl = QL_ALN(ii,jj,5)
          if (neqns>=6) ys01l = QL_ALN(ii,jj,6)
          if (neqns>=7) ys02l = QL_ALN(ii,jj,7)
          if (neqns>=8) ys03l = QL_ALN(ii,jj,8)
          if (neqns>=9) ys04l = QL_ALN(ii,jj,9)
          if (neqns>=10)ys05l = QL_ALN(ii,jj,10)
          if (neqns>=11)ys06l = QL_ALN(ii,jj,11)
          if (neqns>=12)ys07l = QL_ALN(ii,jj,12)

          rhl = DNST_L_ALN(ii,jj)
          hsl = HTPS_L_ALN(ii,jj)

          if (iact_table_thermo==1) then
            rhtl = DRTM_L_ALN(ii,jj)
            rhpl = DRPR_L_ALN(ii,jj)
            hstl = SHCP_L_ALN(ii,jj)*CPREF
            hspl = DHPR_L_ALN(ii,jj)
          else
            rhtl = DRTM_ALN(ii,jl)
            rhpl = DRPR_ALN(ii,jl)
            hstl = SHCP_ALN(ii,jl)*CPREF
            hspl = DHPR_ALN(ii,jl)
          endif

          q01l = tl
          q02l = ul
          q03l = vl
          q04l = wl
          q05l = pl
          if (neqns>=6) q06l = ys01l
          if (neqns>=7) q07l = ys02l
          if (neqns>=8) q08l = ys03l
          if (neqns>=9) q09l = ys04l
          if (neqns>=10)q10l = ys05l
          if (neqns>=11)q11l = ys06l
          if (neqns>=12)q12l = ys07l

          qql = ul*ul+vl*vl+wl*wl

          htl = hsl+0.5d0*qql

          ucl = sxm*ul+sym*vl+szm*wl

          f01l = rhl*(ucl-acm)
          f02l = rhl*(ucl-acm)*ul+sxm*pl
          f03l = rhl*(ucl-acm)*vl+sym*pl
          f04l = rhl*(ucl-acm)*wl+szm*pl
          f05l = rhl*(ucl-acm)*htl+acm*pl
          if (neqns>=6) f06l = rhl*(ucl-acm)*ys01l
          if (neqns>=7) f07l = rhl*(ucl-acm)*ys02l
          if (neqns>=8) f08l = rhl*(ucl-acm)*ys03l
          if (neqns>=9) f09l = rhl*(ucl-acm)*ys04l
          if (neqns>=10)f10l = rhl*(ucl-acm)*ys05l
          if (neqns>=11)f11l = rhl*(ucl-acm)*ys06l
          if (neqns>=12)f12l = rhl*(ucl-acm)*ys07l

          c2l = rhl*hstl/(rhl*rhpl*hstl+rhtl*(1.0d0-rhl*hspl))
          cl  = dsqrt(c2l)

          vpl = dmax1(VPRE_ALN(ii,jl),VPRE_ALN(ii,jr))
          vpl = dmin1(vpl,cl)
          vp2l= vpl*vpl

          rrh  = 1.0d0/rhl
          hsqq = hsl-0.5d0*qql
          alph = 1.0d0-rhl*hspl
          dp   = rhl*rhpl*hstl+rhtl*alph
          rdp  = 1.0d0/dp

          g01l = rdp*(alph-rhpl*hsqq)*f01l-rdp*rhpl*ul*f02l-rdp*rhpl*vl*f03l-rdp*rhpl*wl*f04l+rdp*rhpl*f05l
          g02l =-rrh*ul*f01l+rrh*f02l
          g03l =-rrh*vl*f01l+rrh*f03l
          g04l =-rrh*wl*f01l+rrh*f04l
          g05l = rdp*(rhl*hstl+rhtl*hsqq)*f01l+rdp*rhtl*ul*f02l+rdp*rhtl*vl*f03l+rdp*rhtl*wl*f04l-rdp*rhtl*f05l
          if (neqns>=6) g06l =-rrh*ys01l*f01l+rrh*f06l
          if (neqns>=7) g07l =-rrh*ys02l*f01l+rrh*f07l
          if (neqns>=8) g08l =-rrh*ys03l*f01l+rrh*f08l
          if (neqns>=9) g09l =-rrh*ys04l*f01l+rrh*f09l
          if (neqns>=10)g10l =-rrh*ys05l*f01l+rrh*f10l
          if (neqns>=11)g11l =-rrh*ys06l*f01l+rrh*f11l
          if (neqns>=12)g12l =-rrh*ys07l*f01l+rrh*f12l

          f01l = g01l-alph*(1.0d0-vp2l/c2l)/(rhl*hstl)*g05l
          f02l = g02l
          f03l = g03l
          f04l = g04l
          f05l = vp2l/c2l*g05l
          if (neqns>=6) f06l = g06l
          if (neqns>=7) f07l = g07l
          if (neqns>=8) f08l = g08l
          if (neqns>=9) f09l = g09l
          if (neqns>=10)f10l = g10l
          if (neqns>=11)f11l = g11l
          if (neqns>=12)f12l = g12l

          ggdl= dsqrt((ucl-acm)*(ucl-acm)*(1.0d0-vp2l/c2l)*(1.0d0-vp2l/c2l)+4.0d0*vp2l*sxyzm2)

          eigal = ucl-acm
          eigbl = 0.5d0*((ucl-acm)*(1.0d0+vp2l/c2l)+ggdl)
          eigcl = 0.5d0*((ucl-acm)*(1.0d0+vp2l/c2l)-ggdl)

          eig1l = eigal
          eig2l = (1.0d0-cd1)*eigal+cd1*eigbl
          eig3l = (1.0d0-cd2)*eigal+cd2*eigbl
          eig4l = (1.0d0-cd3)*eigal+cd3*eigbl
          eig5l = eigcl

          tht0 = (1.0d0-rhl*hspl)/(rhl*hstl)
          tht2 = sxyzm2/(rhl*(eigbl-vp2l/c2l*eigal))
          tht3 = sxyzm2/(rhl*(eigcl-vp2l/c2l*eigal))

          thtb = 1.0d0/(tht2-tht3)
          thtb0= tht0*thtb
          thtb2= tht2*thtb
          thtb3= tht3*thtb

          !L
          m11 = 1.0d0
          m12 = 0.0d0
          m13 = 0.0d0
          m14 = 0.0d0
          m15 =-tht0
          !
          m21 = 0.0d0
          m22 = smm
          m23 = cd1*sym-cd2*sxm
          m24 = cd1*szm-cd3*sxm
          m25 = cd1*tht2
          !
          m31 = 0.0d0
          m32 = cd2*sxm-cd1*sym
          m33 = smm
          m34 = cd2*szm-cd3*sym
          m35 = cd2*tht2
          !
          m41 = 0.0d0
          m42 = cd3*sxm-cd1*szm
          m43 = cd3*sym-cd2*szm
          m44 = smm
          m45 = cd3*tht2
          !
          m51 = 0.0d0
          m52 = sxm
          m53 = sym
          m54 = szm
          m55 = tht3

          g01l = m11*f01l + m12*f02l + m13*f03l + m14*f04l + m15*f05l
          g02l = m21*f01l + m22*f02l + m23*f03l + m24*f04l + m25*f05l
          g03l = m31*f01l + m32*f02l + m33*f03l + m34*f04l + m35*f05l
          g04l = m41*f01l + m42*f02l + m43*f03l + m44*f04l + m45*f05l
          g05l = m51*f01l + m52*f02l + m53*f03l + m54*f04l + m55*f05l
          if (neqns>=6) g06l = f06l
          if (neqns>=7) g07l = f07l
          if (neqns>=8) g08l = f08l
          if (neqns>=9) g09l = f09l
          if (neqns>=10)g10l = f10l
          if (neqns>=11)g11l = f11l
          if (neqns>=12)g12l = f12l

          s01l = m11*q01l + m12*q02l + m13*q03l + m14*q04l + m15*q05l
          s02l = m21*q01l + m22*q02l + m23*q03l + m24*q04l + m25*q05l
          s03l = m31*q01l + m32*q02l + m33*q03l + m34*q04l + m35*q05l
          s04l = m41*q01l + m42*q02l + m43*q03l + m44*q04l + m45*q05l
          s05l = m51*q01l + m52*q02l + m53*q03l + m54*q04l + m55*q05l
          if (neqns>=6) s06l = q06l
          if (neqns>=7) s07l = q07l
          if (neqns>=8) s08l = q08l
          if (neqns>=9) s09l = q09l
          if (neqns>=10)s10l = q10l
          if (neqns>=11)s11l = q11l
          if (neqns>=12)s12l = q12l


          !right side
          tr = QR_ALN(ii,jj,1)
          ur = QR_ALN(ii,jj,2)
          vr = QR_ALN(ii,jj,3)
          wr = QR_ALN(ii,jj,4)
          pr = QR_ALN(ii,jj,5)
          if (neqns>=6) ys01r = QR_ALN(ii,jj,6)
          if (neqns>=7) ys02r = QR_ALN(ii,jj,7)
          if (neqns>=8) ys03r = QR_ALN(ii,jj,8)
          if (neqns>=9) ys04r = QR_ALN(ii,jj,9)
          if (neqns>=10)ys05r = QR_ALN(ii,jj,10)
          if (neqns>=11)ys06r = QR_ALN(ii,jj,11)
          if (neqns>=12)ys07r = QR_ALN(ii,jj,12)

          rhr = DNST_R_ALN(ii,jj)
          hsr = HTPS_R_ALN(ii,jj)

          if (iact_table_thermo==1) then
            rhtr = DRTM_R_ALN(ii,jj)
            rhpr = DRPR_R_ALN(ii,jj)
            hstr = SHCP_R_ALN(ii,jj)*CPREF
            hspr = DHPR_R_ALN(ii,jj)
          else
            rhtr = DRTM_ALN(ii,jr)
            rhpr = DRPR_ALN(ii,jr)
            hstr = SHCP_ALN(ii,jr)*CPREF
            hspr = DHPR_ALN(ii,jr)
          endif

          q01r = tr
          q02r = ur
          q03r = vr
          q04r = wr
          q05r = pr
          if (neqns>=6) q06r = ys01r
          if (neqns>=7) q07r = ys02r
          if (neqns>=8) q08r = ys03r
          if (neqns>=9) q09r = ys04r
          if (neqns>=10)q10r = ys05r
          if (neqns>=11)q11r = ys06r
          if (neqns>=12)q12r = ys07r

          qqr = ur*ur+vr*vr+wr*wr

          htr = hsr+0.5d0*qqr

          ucr = sxm*ur+sym*vr+szm*wr

          f01r = rhr*(ucr-acm)
          f02r = rhr*(ucr-acm)*ur+sxm*pr
          f03r = rhr*(ucr-acm)*vr+sym*pr
          f04r = rhr*(ucr-acm)*wr+szm*pr
          f05r = rhr*(ucr-acm)*htr+acm*pr
          if (neqns>=6) f06r = rhr*(ucr-acm)*ys01r
          if (neqns>=7) f07r = rhr*(ucr-acm)*ys02r
          if (neqns>=8) f08r = rhr*(ucr-acm)*ys03r
          if (neqns>=9) f09r = rhr*(ucr-acm)*ys04r
          if (neqns>=10)f10r = rhr*(ucr-acm)*ys05r
          if (neqns>=11)f11r = rhr*(ucr-acm)*ys06r
          if (neqns>=12)f12r = rhr*(ucr-acm)*ys07r

          c2r = rhr*hstr/(rhr*rhpr*hstr+rhtr*(1.0d0-rhr*hspr))
          cr  = dsqrt(c2r)

          vpr = dmax1(VPRE_ALN(ii,jl),VPRE_ALN(ii,jr))
          vpr = dmin1(vpr,cr)
          vp2r= vpr*vpr

          rrh  = 1.0d0/rhr
          hsqq = hsr-0.5d0*qqr
          alph = 1.0d0-rhr*hspr
          dp   = rhr*rhpr*hstr+rhtr*alph
          rdp  = 1.0d0/dp

          g01r = rdp*(alph-rhpr*hsqq)*f01r-rdp*rhpr*ur*f02r-rdp*rhpr*vr*f03r-rdp*rhpr*wr*f04r+rdp*rhpr*f05r
          g02r =-rrh*ur*f01r+rrh*f02r
          g03r =-rrh*vr*f01r+rrh*f03r
          g04r =-rrh*wr*f01r+rrh*f04r
          g05r = rdp*(rhr*hstr+rhtr*hsqq)*f01r+rdp*rhtr*ur*f02r+rdp*rhtr*vr*f03r+rdp*rhtr*wr*f04r-rdp*rhtr*f05r
          if (neqns>=6) g06r =-rrh*ys01r*f01r+rrh*f06r
          if (neqns>=7) g07r =-rrh*ys02r*f01r+rrh*f07r
          if (neqns>=8) g08r =-rrh*ys03r*f01r+rrh*f08r
          if (neqns>=9) g09r =-rrh*ys04r*f01r+rrh*f09r
          if (neqns>=10)g10r =-rrh*ys05r*f01r+rrh*f10r
          if (neqns>=11)g11r =-rrh*ys06r*f01r+rrh*f11r
          if (neqns>=12)g12r =-rrh*ys07r*f01r+rrh*f12r

          f01r = g01r-alph*(1.0d0-vp2r/c2r)/(rhr*hstr)*g05r
          f02r = g02r
          f03r = g03r
          f04r = g04r
          f05r = vp2r/c2r*g05r
          if (neqns>=6) f06r = g06r
          if (neqns>=7) f07r = g07r
          if (neqns>=8) f08r = g08r
          if (neqns>=9) f09r = g09r
          if (neqns>=10)f10r = g10r
          if (neqns>=11)f11r = g11r
          if (neqns>=12)f12r = g12r

          ggdr= dsqrt((ucr-acm)*(ucr-acm)*(1.0d0-vp2r/c2r)*(1.0d0-vp2r/c2r)+4.0d0*vp2r*sxyzm2)

          eigar = ucr-acm
          eigbr = 0.5d0*((ucr-acm)*(1.0d0+vp2r/c2r)+ggdr)
          eigcr = 0.5d0*((ucr-acm)*(1.0d0+vp2r/c2r)-ggdr)

          eig1r = eigar
          eig2r = (1.0d0-cd1)*eigar+cd1*eigbr
          eig3r = (1.0d0-cd2)*eigar+cd2*eigbr
          eig4r = (1.0d0-cd3)*eigar+cd3*eigbr
          eig5r = eigcr

          tht0 = (1.0d0-rhr*hspr)/(rhr*hstr)
          tht2 = sxyzm2/(rhr*(eigbr-vp2r/c2r*eigar))
          tht3 = sxyzm2/(rhr*(eigcr-vp2r/c2r*eigar))

          thtb = 1.0d0/(tht2-tht3)
          thtb0= tht0*thtb
          thtb2= tht2*thtb
          thtb3= tht3*thtb

          !L
          m11 = 1.0d0
          m12 = 0.0d0
          m13 = 0.0d0
          m14 = 0.0d0
          m15 =-tht0
          !
          m21 = 0.0d0
          m22 = smm
          m23 = cd1*sym-cd2*sxm
          m24 = cd1*szm-cd3*sxm
          m25 = cd1*tht2
          !
          m31 = 0.0d0
          m32 = cd2*sxm-cd1*sym
          m33 = smm
          m34 = cd2*szm-cd3*sym
          m35 = cd2*tht2
          !
          m41 = 0.0d0
          m42 = cd3*sxm-cd1*szm
          m43 = cd3*sym-cd2*szm
          m44 = smm
          m45 = cd3*tht2
          !
          m51 = 0.0d0
          m52 = sxm
          m53 = sym
          m54 = szm
          m55 = tht3

          g01r = m11*f01r + m12*f02r + m13*f03r + m14*f04r + m15*f05r
          g02r = m21*f01r + m22*f02r + m23*f03r + m24*f04r + m25*f05r
          g03r = m31*f01r + m32*f02r + m33*f03r + m34*f04r + m35*f05r
          g04r = m41*f01r + m42*f02r + m43*f03r + m44*f04r + m45*f05r
          g05r = m51*f01r + m52*f02r + m53*f03r + m54*f04r + m55*f05r
          if (neqns>=6) g06r = f06r
          if (neqns>=7) g07r = f07r
          if (neqns>=8) g08r = f08r
          if (neqns>=9) g09r = f09r
          if (neqns>=10)g10r = f10r
          if (neqns>=11)g11r = f11r
          if (neqns>=12)g12r = f12r

          s01r = m11*q01r + m12*q02r + m13*q03r + m14*q04r + m15*q05r
          s02r = m21*q01r + m22*q02r + m23*q03r + m24*q04r + m25*q05r
          s03r = m31*q01r + m32*q02r + m33*q03r + m34*q04r + m35*q05r
          s04r = m41*q01r + m42*q02r + m43*q03r + m44*q04r + m45*q05r
          s05r = m51*q01r + m52*q02r + m53*q03r + m54*q04r + m55*q05r
          if (neqns>=6) s06r = q06r
          if (neqns>=7) s07r = q07r
          if (neqns>=8) s08r = q08r
          if (neqns>=9) s09r = q09r
          if (neqns>=10)s10r = q10r
          if (neqns>=11)s11r = q11r
          if (neqns>=12)s12r = q12r

          eig1max = dmax1(dabs(eig1l),dabs(eig1r))
          eig2max = dmax1(dabs(eig2l),dabs(eig2r))
          eig3max = dmax1(dabs(eig3l),dabs(eig3r))
          eig4max = dmax1(dabs(eig4l),dabs(eig4r))
          eig5max = dmax1(dabs(eig5l),dabs(eig5r))

          llf01l = 0.5d0*(g01l+eig1max*s01l)
          llf02l = 0.5d0*(g02l+eig2max*s02l)
          llf03l = 0.5d0*(g03l+eig3max*s03l)
          llf04l = 0.5d0*(g04l+eig4max*s04l)
          llf05l = 0.5d0*(g05l+eig5max*s05l)
          if (neqns>=6) llf06l = 0.5d0*(g06l+eig1max*s06l)
          if (neqns>=7) llf07l = 0.5d0*(g07l+eig1max*s07l)
          if (neqns>=8) llf08l = 0.5d0*(g08l+eig1max*s08l)
          if (neqns>=9) llf09l = 0.5d0*(g09l+eig1max*s09l)
          if (neqns>=10)llf10l = 0.5d0*(g10l+eig1max*s10l)
          if (neqns>=11)llf11l = 0.5d0*(g11l+eig1max*s11l)
          if (neqns>=12)llf12l = 0.5d0*(g12l+eig1max*s12l)

          llf01r = 0.5d0*(g01r-eig1max*s01r)
          llf02r = 0.5d0*(g02r-eig2max*s02r)
          llf03r = 0.5d0*(g03r-eig3max*s03r)
          llf04r = 0.5d0*(g04r-eig4max*s04r)
          llf05r = 0.5d0*(g05r-eig5max*s05r)
          if (neqns>=6) llf06r = 0.5d0*(g06r-eig1max*s06r)
          if (neqns>=7) llf07r = 0.5d0*(g07r-eig1max*s07r)
          if (neqns>=8) llf08r = 0.5d0*(g08r-eig1max*s08r)
          if (neqns>=9) llf09r = 0.5d0*(g09r-eig1max*s09r)
          if (neqns>=10)llf10r = 0.5d0*(g10r-eig1max*s10r)
          if (neqns>=11)llf11r = 0.5d0*(g11r-eig1max*s11r)
          if (neqns>=12)llf12r = 0.5d0*(g12r-eig1max*s12r)

          sig1 = dsign(1.0d0,eig1l)*dsign(1.0d0,eig1r)
          sig2 = dsign(1.0d0,eig2l)*dsign(1.0d0,eig2r)
          sig3 = dsign(1.0d0,eig3l)*dsign(1.0d0,eig3r)
          sig4 = dsign(1.0d0,eig4l)*dsign(1.0d0,eig4r)
          sig5 = dsign(1.0d0,eig5l)*dsign(1.0d0,eig5r)

          sig1l = 0.5d0*(1.0d0+dsign(1.0d0,eig1l))
          sig2l = 0.5d0*(1.0d0+dsign(1.0d0,eig2l))
          sig3l = 0.5d0*(1.0d0+dsign(1.0d0,eig3l))
          sig4l = 0.5d0*(1.0d0+dsign(1.0d0,eig4l))
          sig5l = 0.5d0*(1.0d0+dsign(1.0d0,eig5l))

          sig1r = 0.5d0*(1.0d0-dsign(1.0d0,eig1r))
          sig2r = 0.5d0*(1.0d0-dsign(1.0d0,eig2r))
          sig3r = 0.5d0*(1.0d0-dsign(1.0d0,eig3r))
          sig4r = 0.5d0*(1.0d0-dsign(1.0d0,eig4r))
          sig5r = 0.5d0*(1.0d0-dsign(1.0d0,eig5r))

          xai01l = 0.5d0*((1.0d0+sig1)*sig1l*g01l+(1.0d0-sig1)*llf01l)
          xai02l = 0.5d0*((1.0d0+sig2)*sig2l*g02l+(1.0d0-sig2)*llf02l)
          xai03l = 0.5d0*((1.0d0+sig3)*sig3l*g03l+(1.0d0-sig3)*llf03l)
          xai04l = 0.5d0*((1.0d0+sig4)*sig4l*g04l+(1.0d0-sig4)*llf04l)
          xai05l = 0.5d0*((1.0d0+sig5)*sig5l*g05l+(1.0d0-sig5)*llf05l)
          if (neqns>=6) xai06l = 0.5d0*((1.0d0+sig1)*sig1l*g06l+(1.0d0-sig1)*llf06l)
          if (neqns>=7) xai07l = 0.5d0*((1.0d0+sig1)*sig1l*g07l+(1.0d0-sig1)*llf07l)
          if (neqns>=8) xai08l = 0.5d0*((1.0d0+sig1)*sig1l*g08l+(1.0d0-sig1)*llf08l)
          if (neqns>=9) xai09l = 0.5d0*((1.0d0+sig1)*sig1l*g09l+(1.0d0-sig1)*llf09l)
          if (neqns>=10)xai10l = 0.5d0*((1.0d0+sig1)*sig1l*g10l+(1.0d0-sig1)*llf10l)
          if (neqns>=11)xai11l = 0.5d0*((1.0d0+sig1)*sig1l*g11l+(1.0d0-sig1)*llf11l)
          if (neqns>=12)xai12l = 0.5d0*((1.0d0+sig1)*sig1l*g12l+(1.0d0-sig1)*llf12l)

          xai01r = 0.5d0*((1.0d0+sig1)*sig1r*g01r+(1.0d0-sig1)*llf01r)
          xai02r = 0.5d0*((1.0d0+sig2)*sig2r*g02r+(1.0d0-sig2)*llf02r)
          xai03r = 0.5d0*((1.0d0+sig3)*sig3r*g03r+(1.0d0-sig3)*llf03r)
          xai04r = 0.5d0*((1.0d0+sig4)*sig4r*g04r+(1.0d0-sig4)*llf04r)
          xai05r = 0.5d0*((1.0d0+sig5)*sig5r*g05r+(1.0d0-sig5)*llf05r)
          if (neqns>=6) xai06r = 0.5d0*((1.0d0+sig1)*sig1r*g06r+(1.0d0-sig1)*llf06r)
          if (neqns>=7) xai07r = 0.5d0*((1.0d0+sig1)*sig1r*g07r+(1.0d0-sig1)*llf07r)
          if (neqns>=8) xai08r = 0.5d0*((1.0d0+sig1)*sig1r*g08r+(1.0d0-sig1)*llf08r)
          if (neqns>=9) xai09r = 0.5d0*((1.0d0+sig1)*sig1r*g09r+(1.0d0-sig1)*llf09r)
          if (neqns>=10)xai10r = 0.5d0*((1.0d0+sig1)*sig1r*g10r+(1.0d0-sig1)*llf10r)
          if (neqns>=11)xai11r = 0.5d0*((1.0d0+sig1)*sig1r*g11r+(1.0d0-sig1)*llf11r)
          if (neqns>=12)xai12r = 0.5d0*((1.0d0+sig1)*sig1r*g12r+(1.0d0-sig1)*llf12r)


          !Left side
          tht0 = (1.0d0-rhl*hspl)/(rhl*hstl)
          tht2 = sxyzm2/(rhl*(eigbl-vp2l/c2l*eigal))
          tht3 = sxyzm2/(rhl*(eigcl-vp2l/c2l*eigal))

          thtb = 1.0d0/(tht2-tht3)
          thtb0= tht0*thtb
          thtb2= tht2*thtb
          thtb3= tht3*thtb

          !R
          m11 = 1.0d0
          m12 = thtb0*cd1
          m13 = thtb0*cd2
          m14 = thtb0*cd3
          m15 =-thtb0
          !
          m21 = 0.0d0
          m22 =((1.0d0-thtb3)*cd1-1.0d0)*sxm*sxm/(sxyzm2*smm)+(1.0d0-cd1)/smm
          m23 =((1.0d0-thtb3)*cd2-1.0d0)*sxm*sym/(sxyzm2*smm)
          m24 =((1.0d0-thtb3)*cd3-1.0d0)*sxm*szm/(sxyzm2*smm)
          m25 = thtb2*sxm/sxyzm2
          !
          m31 = 0.0d0
          m32 =((1.0d0-thtb3)*cd1-1.0d0)*sym*sxm/(sxyzm2*smm)
          m33 =((1.0d0-thtb3)*cd2-1.0d0)*sym*sym/(sxyzm2*smm)+(1.0d0-cd2)/smm
          m34 =((1.0d0-thtb3)*cd3-1.0d0)*sym*szm/(sxyzm2*smm)
          m35 = thtb2*sym/sxyzm2
          !
          m41 = 0.0d0
          m42 =((1.0d0-thtb3)*cd1-1.0d0)*szm*sxm/(sxyzm2*smm)
          m43 =((1.0d0-thtb3)*cd2-1.0d0)*szm*sym/(sxyzm2*smm)
          m44 =((1.0d0-thtb3)*cd3-1.0d0)*szm*szm/(sxyzm2*smm)+(1.0d0-cd3)/smm
          m45 = thtb2*szm/sxyzm2
          !
          m51 = 0.0d0
          m52 = thtb*cd1
          m53 = thtb*cd2
          m54 = thtb*cd3
          m55 =-thtb

          f01l = m11*xai01l + m12*xai02l + m13*xai03l + m14*xai04l + m15*xai05l
          f02l = m21*xai01l + m22*xai02l + m23*xai03l + m24*xai04l + m25*xai05l
          f03l = m31*xai01l + m32*xai02l + m33*xai03l + m34*xai04l + m35*xai05l
          f04l = m41*xai01l + m42*xai02l + m43*xai03l + m44*xai04l + m45*xai05l
          f05l = m51*xai01l + m52*xai02l + m53*xai03l + m54*xai04l + m55*xai05l
          if (neqns>=6) f06l = xai06l
          if (neqns>=7) f07l = xai07l
          if (neqns>=8) f08l = xai08l
          if (neqns>=9) f09l = xai09l
          if (neqns>=10)f10l = xai10l
          if (neqns>=11)f11l = xai11l
          if (neqns>=12)f12l = xai12l

          alph = 1.0d0-rhl*hspl

          g01l = f01l+alph*(c2l/vp2l-1.0d0)/(rhl*hstl)*f05l
          g02l = f02l
          g03l = f03l
          g04l = f04l
          g05l = c2l/vp2l*f05l
          if (neqns>=6) g06l = f06l
          if (neqns>=7) g07l = f07l
          if (neqns>=8) g08l = f08l
          if (neqns>=9) g09l = f09l
          if (neqns>=10)g10l = f10l
          if (neqns>=11)g11l = f11l
          if (neqns>=12)g12l = f12l

          f01l =    rhtl*g01l+rhpl*g05l
          f02l = ul*rhtl*g01l+rhl*g02l+ul*rhpl*g05l
          f03l = vl*rhtl*g01l+rhl*g03l+vl*rhpl*g05l
          f04l = wl*rhtl*g01l+rhl*g04l+wl*rhpl*g05l
          f05l = (rhtl*htl+rhl*hstl)*g01l+ul*rhl*g02l+vl*rhl*g03l+wl*rhl*g04l+(htl*rhpl-alph)*g05l
          if (neqns>=6) f06l = ys01l*rhtl*g01l+ys01l*rhpl*g05l+rhl*g06l
          if (neqns>=7) f07l = ys02l*rhtl*g01l+ys02l*rhpl*g05l+rhl*g07l
          if (neqns>=8) f08l = ys03l*rhtl*g01l+ys03l*rhpl*g05l+rhl*g08l
          if (neqns>=9) f09l = ys04l*rhtl*g01l+ys04l*rhpl*g05l+rhl*g09l
          if (neqns>=10)f10l = ys05l*rhtl*g01l+ys05l*rhpl*g05l+rhl*g10l
          if (neqns>=11)f11l = ys06l*rhtl*g01l+ys06l*rhpl*g05l+rhl*g11l
          if (neqns>=12)f12l = ys07l*rhtl*g01l+ys07l*rhpl*g05l+rhl*g12l


          !Right side
          tht0 = (1.0d0-rhr*hspr)/(rhr*hstr)
          tht2 = sxyzm2/(rhr*(eigbr-vp2r/c2r*eigar))
          tht3 = sxyzm2/(rhr*(eigcr-vp2r/c2r*eigar))

          thtb = 1.0d0/(tht2-tht3)
          thtb0= tht0*thtb
          thtb2= tht2*thtb
          thtb3= tht3*thtb

          !R
          m11 = 1.0d0
          m12 = thtb0*cd1
          m13 = thtb0*cd2
          m14 = thtb0*cd3
          m15 =-thtb0
          !
          m21 = 0.0d0
          m22 =((1.0d0-thtb3)*cd1-1.0d0)*sxm*sxm/(sxyzm2*smm)+(1.0d0-cd1)/smm
          m23 =((1.0d0-thtb3)*cd2-1.0d0)*sxm*sym/(sxyzm2*smm)
          m24 =((1.0d0-thtb3)*cd3-1.0d0)*sxm*szm/(sxyzm2*smm)
          m25 = thtb2*sxm/sxyzm2
          !
          m31 = 0.0d0
          m32 =((1.0d0-thtb3)*cd1-1.0d0)*sym*sxm/(sxyzm2*smm)
          m33 =((1.0d0-thtb3)*cd2-1.0d0)*sym*sym/(sxyzm2*smm)+(1.0d0-cd2)/smm
          m34 =((1.0d0-thtb3)*cd3-1.0d0)*sym*szm/(sxyzm2*smm)
          m35 = thtb2*sym/sxyzm2
          !
          m41 = 0.0d0
          m42 =((1.0d0-thtb3)*cd1-1.0d0)*szm*sxm/(sxyzm2*smm)
          m43 =((1.0d0-thtb3)*cd2-1.0d0)*szm*sym/(sxyzm2*smm)
          m44 =((1.0d0-thtb3)*cd3-1.0d0)*szm*szm/(sxyzm2*smm)+(1.0d0-cd3)/smm
          m45 = thtb2*szm/sxyzm2
          !
          m51 = 0.0d0
          m52 = thtb*cd1
          m53 = thtb*cd2
          m54 = thtb*cd3
          m55 =-thtb

          f01r = m11*xai01r + m12*xai02r + m13*xai03r + m14*xai04r + m15*xai05r
          f02r = m21*xai01r + m22*xai02r + m23*xai03r + m24*xai04r + m25*xai05r
          f03r = m31*xai01r + m32*xai02r + m33*xai03r + m34*xai04r + m35*xai05r
          f04r = m41*xai01r + m42*xai02r + m43*xai03r + m44*xai04r + m45*xai05r
          f05r = m51*xai01r + m52*xai02r + m53*xai03r + m54*xai04r + m55*xai05r
          if (neqns>=6) f06r = xai06r
          if (neqns>=7) f07r = xai07r
          if (neqns>=8) f08r = xai08r
          if (neqns>=9) f09r = xai09r
          if (neqns>=10)f10r = xai10r
          if (neqns>=11)f11r = xai11r
          if (neqns>=12)f12r = xai12r

          alph = 1.0d0-rhr*hspr

          g01r = f01r+alph*(c2r/vp2r-1.0d0)/(rhr*hstr)*f05r
          g02r = f02r
          g03r = f03r
          g04r = f04r
          g05r = c2r/vp2r*f05r
          if (neqns>=6) g06r = f06r
          if (neqns>=7) g07r = f07r
          if (neqns>=8) g08r = f08r
          if (neqns>=9) g09r = f09r
          if (neqns>=10)g10r = f10r
          if (neqns>=11)g11r = f11r
          if (neqns>=12)g12r = f12r

          f01r =    rhtr*g01r+rhpr*g05r
          f02r = ur*rhtr*g01r+rhr*g02r+ur*rhpr*g05r
          f03r = vr*rhtr*g01r+rhr*g03r+vr*rhpr*g05r
          f04r = wr*rhtr*g01r+rhr*g04r+wr*rhpr*g05r
          f05r = (rhtr*htr+rhr*hstr)*g01r+ur*rhr*g02r+vr*rhr*g03r+wr*rhr*g04r+(htr*rhpr-alph)*g05r
          if (neqns>=6) f06r = ys01r*rhtr*g01r+ys01r*rhpr*g05r+rhr*g06r
          if (neqns>=7) f07r = ys02r*rhtr*g01r+ys02r*rhpr*g05r+rhr*g07r
          if (neqns>=8) f08r = ys03r*rhtr*g01r+ys03r*rhpr*g05r+rhr*g08r
          if (neqns>=9) f09r = ys04r*rhtr*g01r+ys04r*rhpr*g05r+rhr*g09r
          if (neqns>=10)f10r = ys05r*rhtr*g01r+ys05r*rhpr*g05r+rhr*g10r
          if (neqns>=11)f11r = ys06r*rhtr*g01r+ys06r*rhpr*g05r+rhr*g11r
          if (neqns>=12)f12r = ys07r*rhtr*g01r+ys07r*rhpr*g05r+rhr*g12r


          FLUX_ALN(ii,jj,1) = f01l+f01r
          FLUX_ALN(ii,jj,2) = f02l+f02r
          FLUX_ALN(ii,jj,3) = f03l+f03r
          FLUX_ALN(ii,jj,4) = f04l+f04r
          FLUX_ALN(ii,jj,5) = f05l+f05r
          if (neqns>=6) FLUX_ALN(ii,jj,6) = f06l+f06r
          if (neqns>=7) FLUX_ALN(ii,jj,7) = f07l+f07r
          if (neqns>=8) FLUX_ALN(ii,jj,8) = f08l+f08r
          if (neqns>=9) FLUX_ALN(ii,jj,9) = f09l+f09r
          if (neqns>=10)FLUX_ALN(ii,jj,10)= f10l+f10r
          if (neqns>=11)FLUX_ALN(ii,jj,11)= f11l+f11r
          if (neqns>=12)FLUX_ALN(ii,jj,12)= f12l+f12r
        enddo
      enddo

      do jj=1,jjmax
        do ii=1,iimax
          DFLX_ALN(ii,jj,1) =-(FLUX_ALN(ii,jj+1,1)-FLUX_ALN(ii,jj,1))
          DFLX_ALN(ii,jj,2) =-(FLUX_ALN(ii,jj+1,2)-FLUX_ALN(ii,jj,2))
          DFLX_ALN(ii,jj,3) =-(FLUX_ALN(ii,jj+1,3)-FLUX_ALN(ii,jj,3))
          DFLX_ALN(ii,jj,4) =-(FLUX_ALN(ii,jj+1,4)-FLUX_ALN(ii,jj,4))
          DFLX_ALN(ii,jj,5) =-(FLUX_ALN(ii,jj+1,5)-FLUX_ALN(ii,jj,5))
          if (neqns>=6) DFLX_ALN(ii,jj,6) =-(FLUX_ALN(ii,jj+1, 6)-FLUX_ALN(ii,jj, 6))
          if (neqns>=7) DFLX_ALN(ii,jj,7) =-(FLUX_ALN(ii,jj+1, 7)-FLUX_ALN(ii,jj, 7))
          if (neqns>=8) DFLX_ALN(ii,jj,8) =-(FLUX_ALN(ii,jj+1, 8)-FLUX_ALN(ii,jj, 8))
          if (neqns>=9) DFLX_ALN(ii,jj,9) =-(FLUX_ALN(ii,jj+1, 9)-FLUX_ALN(ii,jj, 9))
          if (neqns>=10)DFLX_ALN(ii,jj,10)=-(FLUX_ALN(ii,jj+1,10)-FLUX_ALN(ii,jj,10))
          if (neqns>=11)DFLX_ALN(ii,jj,11)=-(FLUX_ALN(ii,jj+1,11)-FLUX_ALN(ii,jj,11))
          if (neqns>=12)DFLX_ALN(ii,jj,12)=-(FLUX_ALN(ii,jj+1,12)-FLUX_ALN(ii,jj,12))
        enddo
      enddo

    end subroutine calc_fds

  end subroutine flux_riemann_fds_real_gas_extension

end module mdl_flux_riemann_fds_real_gas_extension
