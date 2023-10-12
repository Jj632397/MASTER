module mdl_flux_riemann_fvs_real_gas_extension
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
!Steger and Warming's flux vector splitting scheme with MUSCL or 5th-order WENO.
!Real gas extended version.
!-------------------------------------------------------------------------------

  integer, private, parameter :: isc_reconstruction = 1 !<1>:weno5 <else>:MUSCL

  integer, private, parameter :: iorder_muscl = 3

  integer, private, parameter :: iact_table_thermo = 0

  !Additional numerical dissipation on temperature
  real(8), private, parameter :: coef_thermo_diss = 0.5d0

contains

  subroutine flux_riemann_fvs_real_gas_extension(BLK,RHS)
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

      call calc_fvs(1.0d0,0.0d0,0.0d0)

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

      call calc_fvs(0.0d0,1.0d0,0.0d0)

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

      call calc_fvs(0.0d0,0.0d0,1.0d0)

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
    subroutine calc_fvs(cd1,cd2,cd3)
      real(8), intent(in) :: cd1,cd2,cd3

      integer :: ii,jj,jl,jr,l
      real(8) :: sxm,sym,szm,sxyzm2,sxyzm,smm
      real(8) :: rhl,ul,vl,wl,hsl,htl,qql,el,rhtl,rhpl,hstl,hspl,tl,pl,c2l,cl,ucl,vpl,vp2l,ggdl
      real(8) :: ys1l,ys2l,ys3l,ys4l,ys5l,ys6l,ys7l
      real(8) :: q1l,q2l,q3l,q4l,q5l,q6l,q7l,q8l,q9l,q10l,q11l,q12l
      real(8) :: f1l,f2l,f3l,f4l,f5l,f6l,f7l,f8l,f9l,f10l,f11l,f12l
      real(8) :: eigal,eigbl,eigcl
      real(8) :: eig1l,eig2l,eig3l,eig4l,eig5l
      real(8) :: rhr,ur,vr,wr,hsr,htr,qqr,er,rhtr,rhpr,hstr,hspr,tr,pr,c2r,cr,ucr,vpr,vp2r,ggdr
      real(8) :: ys1r,ys2r,ys3r,ys4r,ys5r,ys6r,ys7r
      real(8) :: q1r,q2r,q3r,q4r,q5r,q6r,q7r,q8r,q9r,q10r,q11r,q12r
      real(8) :: f1r,f2r,f3r,f4r,f5r,f6r,f7r,f8r,f9r,f10r,f11r,f12r
      real(8) :: eigar,eigbr,eigcr
      real(8) :: eig1r,eig2r,eig3r,eig4r,eig5r
      real(8) :: rhm,um,vm,wm,htm,qqm,em,rhtm,rhpm,hstm,hspm,pm,hm,c2m,cm,ucm,vpm,vp2m,ggdm
      real(8) :: eigam,eigbm,eigcm
      real(8) :: eig1m,eig2m,eig3m,eig4m,eig5m
      real(8) :: eig01max,eig02max,eig03max,eig04max,eig05max,eigmax
      real(8) :: eig06max,eig07max,eig08max,eig09max,eig10max,eig11max,eig12max
      real(8) :: dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,dq9,dq10,dq11,dq12
      real(8) :: ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,ds9,ds10,ds11,ds12
      real(8) :: tht0,tht2,tht3,thtb,thtb0,thtb2,thtb3

      real(8) :: m11,m12,m13,m14,m15
      real(8) :: m21,m22,m23,m24,m25
      real(8) :: m31,m32,m33,m34,m35
      real(8) :: m41,m42,m43,m44,m45
      real(8) :: m51,m52,m53,m54,m55

      real(8) :: rhb1,rhb2,rhb3,rhb4,rhb5,rhb6,rhb7,rhb8,rhb9,rhb10,rhb11,rhb12
      real(8) :: rhc1,rhc2,rhc3,rhc4,rhc5,rhc6,rhc7,rhc8,rhc9,rhc10,rhc11,rhc12

      real(8) :: rh,u,v,w,hs,rht,rhp,hst,hsp,c2,c,vp,vp2,qq,ht,alph
      real(8) :: ys1,ys2,ys3,ys4,ys5,ys6,ys7

      real(8) :: pi
      real(8) :: nx,ny,nz,ox,oy,oz
      real(8) :: am,bm,gm,acm,rxm,rym,rzm

      real(8), allocatable, dimension(:,:) :: DTMP
      real(8), allocatable, dimension(:,:) :: THTA

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

      allocate(DTMP(iimax,jjmax+1))
      allocate(THTA(iimax,jjmax  ))

      do jj=1,jjmax+1
        do ii=1,iimax
          DTMP(ii,jj) = dabs(QR_ALN(ii,jj,1)-QL_ALN(ii,jj,1))
        enddo
      enddo

      do jj=1,jjmax
        jl = jj
        jr = jj+1
        do ii=1,iimax
          THTA(ii,jj) = dabs(DTMP(ii,jr)-DTMP(ii,jl))/(DTMP(ii,jr)+DTMP(ii,jl)+1.0d-7)
        enddo
      enddo

      do jj=1,jjmax+1
        jl = max0(jj-1,1    )
        jr = min0(jj  ,jjmax)
        do ii=1,iimax
          DTMP(ii,jj) = coef_thermo_diss*dmax1(THTA(ii,jl),THTA(ii,jr))
        enddo
      enddo

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
          if (neqns>=6) ys1l= QL_ALN(ii,jj,6)
          if (neqns>=7) ys2l= QL_ALN(ii,jj,7)
          if (neqns>=8) ys3l= QL_ALN(ii,jj,8)
          if (neqns>=9) ys4l= QL_ALN(ii,jj,9)
          if (neqns>=10)ys5l= QL_ALN(ii,jj,10)
          if (neqns>=11)ys6l= QL_ALN(ii,jj,11)
          if (neqns>=12)ys7l= QL_ALN(ii,jj,12)

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

          qql = ul*ul+vl*vl+wl*wl

          htl = hsl+0.5d0*qql

          el = rhl*htl-pl

          c2l = rhl*hstl/(rhl*rhpl*hstl+rhtl*(1.0d0-rhl*hspl))
          cl  = dsqrt(c2l)

          vpl = dmin1(VPRE_ALN(ii,jl),VPRE_ALN(ii,jr))
          vpl = dmin1(vpl,cl)
          vp2l= vpl*vpl

          ucl = sxm*ul+sym*vl+szm*wl

          ggdl= dsqrt((ucl-acm)*(ucl-acm)*(1.0d0-vp2l/c2l)*(1.0d0-vp2l/c2l)+4.0d0*vp2l*sxyzm2)

          q1l = tl
          q2l = ul
          q3l = vl
          q4l = wl
          q5l = pl
          if (neqns>=6) q6l = ys1l
          if (neqns>=7) q7l = ys2l
          if (neqns>=8) q8l = ys3l
          if (neqns>=9) q9l = ys4l
          if (neqns>=10)q10l= ys5l
          if (neqns>=11)q11l= ys6l
          if (neqns>=12)q12l= ys7l

          f1l = rhl*(ucl-acm)
          f2l = rhl*(ucl-acm)*ul+sxm*pl
          f3l = rhl*(ucl-acm)*vl+sym*pl
          f4l = rhl*(ucl-acm)*wl+szm*pl
          f5l = rhl*(ucl-acm)*htl+acm*pl
          if (neqns>=6) f6l = rhl*(ucl-acm)*ys1l
          if (neqns>=7) f7l = rhl*(ucl-acm)*ys2l
          if (neqns>=8) f8l = rhl*(ucl-acm)*ys3l
          if (neqns>=9) f9l = rhl*(ucl-acm)*ys4l
          if (neqns>=10)f10l= rhl*(ucl-acm)*ys5l
          if (neqns>=11)f11l= rhl*(ucl-acm)*ys6l
          if (neqns>=12)f12l= rhl*(ucl-acm)*ys7l

          eigal = ucl-acm
          eigbl = 0.5d0*((ucl-acm)*(1.0d0+vp2l/c2l)+ggdl)
          eigcl = 0.5d0*((ucl-acm)*(1.0d0+vp2l/c2l)-ggdl)

          eig1l = eigal
          eig2l = (1.0d0-cd1)*eigal+cd1*eigbl
          eig3l = (1.0d0-cd2)*eigal+cd2*eigbl
          eig4l = (1.0d0-cd3)*eigal+cd3*eigbl
          eig5l = eigcl

          !right side
          tr = QR_ALN(ii,jj,1)
          ur = QR_ALN(ii,jj,2)
          vr = QR_ALN(ii,jj,3)
          wr = QR_ALN(ii,jj,4)
          pr = QR_ALN(ii,jj,5)
          if (neqns>=6) ys1r= QR_ALN(ii,jj,6)
          if (neqns>=7) ys2r= QR_ALN(ii,jj,7)
          if (neqns>=8) ys3r= QR_ALN(ii,jj,8)
          if (neqns>=9) ys4r= QR_ALN(ii,jj,9)
          if (neqns>=10)ys5r= QR_ALN(ii,jj,10)
          if (neqns>=11)ys6r= QR_ALN(ii,jj,11)
          if (neqns>=12)ys7r= QR_ALN(ii,jj,12)

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

          qqr = ur*ur+vr*vr+wr*wr

          htr = hsr+0.5d0*qqr

          er = rhr*htr-pr

          c2r = rhr*hstr/(rhr*rhpr*hstr+rhtr*(1.0d0-rhr*hspr))
          cr  = dsqrt(c2r)

          vpr = dmin1(VPRE_ALN(ii,jl),VPRE_ALN(ii,jr))
          vpr = dmin1(vpr,cr)
          vp2r= vpr*vpr

          ucr = sxm*ur+sym*vr+szm*wr

          ggdr= dsqrt((ucr-acm)*(ucr-acm)*(1.0d0-vp2r/c2r)*(1.0d0-vp2r/c2r)+4.0d0*vp2r*sxyzm2)

          q1r = tr
          q2r = ur
          q3r = vr
          q4r = wr
          q5r = pr
          if (neqns>=6) q6r = ys1r
          if (neqns>=7) q7r = ys2r
          if (neqns>=8) q8r = ys3r
          if (neqns>=9) q9r = ys4r
          if (neqns>=10)q10r= ys5r
          if (neqns>=11)q11r= ys6r
          if (neqns>=12)q12r= ys7r

          f1r = rhr*(ucr-acm)
          f2r = rhr*(ucr-acm)*ur+sxm*pr
          f3r = rhr*(ucr-acm)*vr+sym*pr
          f4r = rhr*(ucr-acm)*wr+szm*pr
          f5r = rhr*(ucr-acm)*htr+acm*pr
          if (neqns>=6) f6r = rhr*(ucr-acm)*ys1r
          if (neqns>=7) f7r = rhr*(ucr-acm)*ys2r
          if (neqns>=8) f8r = rhr*(ucr-acm)*ys3r
          if (neqns>=9) f9r = rhr*(ucr-acm)*ys4r
          if (neqns>=10)f10r= rhr*(ucr-acm)*ys5r
          if (neqns>=11)f11r= rhr*(ucr-acm)*ys6r
          if (neqns>=12)f12r= rhr*(ucr-acm)*ys7r

          eigar = ucr-acm
          eigbr = 0.5d0*((ucr-acm)*(1.0d0+vp2r/c2r)+ggdr)
          eigcr = 0.5d0*((ucr-acm)*(1.0d0+vp2r/c2r)-ggdr)

          eig1r = eigar
          eig2r = (1.0d0-cd1)*eigar+cd1*eigbr
          eig3r = (1.0d0-cd2)*eigar+cd2*eigbr
          eig4r = (1.0d0-cd3)*eigar+cd3*eigbr
          eig5r = eigcr

          !midpoint
          rhm = 0.5d0*(rhl+rhr)
          um  = 0.5d0*( ul+ ur)
          vm  = 0.5d0*( vl+ vr)
          wm  = 0.5d0*( wl+ wr)
          htm = 0.5d0*(htl+htr)

          rhtm = 0.5d0*(rhtl+rhtr)
          rhpm = 0.5d0*(rhpl+rhpr)
          hstm = 0.5d0*(hstl+hstr)
          hspm = 0.5d0*(hspl+hspr)

          c2m = rhm*hstm/(rhm*rhpm*hstm+rhtm*(1.0d0-rhm*hspm))
          cm  = dsqrt(c2m)

          vpm  = dmin1(VPRE_ALN(ii,jl),VPRE_ALN(ii,jr))
          vpm = dmin1(vpm,cm)
          vp2m= vpm*vpm

          ucm = sxm*um+sym*vm+szm*wm

          ggdm= dsqrt((ucm-acm)*(ucm-acm)*(1.0d0-vp2m/c2m)*(1.0d0-vp2m/c2m)+4.0d0*vp2m*sxyzm2)

          eigam = ucm-acm
          eigbm = 0.5d0*((ucm-acm)*(1.0d0+vp2m/c2m)+ggdm)
          eigcm = 0.5d0*((ucm-acm)*(1.0d0+vp2m/c2m)-ggdm)

          eig1m = eigam
          eig2m = (1.0d0-cd1)*eigam+cd1*eigbm
          eig3m = (1.0d0-cd2)*eigam+cd2*eigbm
          eig4m = (1.0d0-cd3)*eigam+cd3*eigbm
          eig5m = eigcm

          !maximum eigenvalues
          eig01max = dmax1(dabs(eig1l),dabs(eig1m),dabs(eig1r))
          eig02max = dmax1(dabs(eig2l),dabs(eig2m),dabs(eig2r))
          eig03max = dmax1(dabs(eig3l),dabs(eig3m),dabs(eig3r))
          eig04max = dmax1(dabs(eig4l),dabs(eig4m),dabs(eig4r))
          eig05max = dmax1(dabs(eig5l),dabs(eig5m),dabs(eig5r))
          !
          if (neqns>=6) eig06max = eig01max
          if (neqns>=7) eig07max = eig01max
          if (neqns>=8) eig08max = eig01max
          if (neqns>=9) eig09max = eig01max
          if (neqns>=10)eig10max = eig01max
          if (neqns>=11)eig11max = eig01max
          if (neqns>=12)eig12max = eig01max

          eigmax = dmax1(eig01max,eig02max,eig03max,eig04max,eig05max)

          eig01max = (1.0d0-DTMP(ii,jj))*eig01max + DTMP(ii,jj)*eigmax

          dq1 = q1r-q1l
          dq2 = q2r-q2l
          dq3 = q3r-q3l
          dq4 = q4r-q4l
          dq5 = q5r-q5l
          if (neqns>=6) dq6 = q6r-q6l
          if (neqns>=7) dq7 = q7r-q7l
          if (neqns>=8) dq8 = q8r-q8l
          if (neqns>=9) dq9 = q9r-q9l
          if (neqns>=10)dq10= q10r-q10l
          if (neqns>=11)dq11= q11r-q11l
          if (neqns>=12)dq12= q12r-q12l

          tht0 = (1.0d0-rhm*hspm)/(rhm*hstm)
          tht2 = sxyzm2/(rhm*(eigbm-vp2m/c2m*eigam))
          tht3 = sxyzm2/(rhm*(eigcm-vp2m/c2m*eigam))

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

          ds1 = m11*dq1 + m12*dq2 + m13*dq3 + m14*dq4 + m15*dq5
          ds2 = m21*dq1 + m22*dq2 + m23*dq3 + m24*dq4 + m25*dq5
          ds3 = m31*dq1 + m32*dq2 + m33*dq3 + m34*dq4 + m35*dq5
          ds4 = m41*dq1 + m42*dq2 + m43*dq3 + m44*dq4 + m45*dq5
          ds5 = m51*dq1 + m52*dq2 + m53*dq3 + m54*dq4 + m55*dq5
          if (neqns>=6) ds6 = dq6
          if (neqns>=7) ds7 = dq7
          if (neqns>=8) ds8 = dq8
          if (neqns>=9) ds9 = dq9
          if (neqns>=10)ds10= dq10
          if (neqns>=11)ds11= dq11
          if (neqns>=12)ds12= dq12

          ds1 = eig01max*ds1
          ds2 = eig02max*ds2
          ds3 = eig03max*ds3
          ds4 = eig04max*ds4
          ds5 = eig05max*ds5
          if (neqns>=6) ds6 = eig06max*ds6
          if (neqns>=7) ds7 = eig07max*ds7
          if (neqns>=8) ds8 = eig08max*ds8
          if (neqns>=9) ds9 = eig09max*ds9
          if (neqns>=10)ds10= eig10max*ds10
          if (neqns>=11)ds11= eig11max*ds11
          if (neqns>=12)ds12= eig12max*ds12

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

          dq1 = m11*ds1 + m12*ds2 + m13*ds3 + m14*ds4 + m15*ds5
          dq2 = m21*ds1 + m22*ds2 + m23*ds3 + m24*ds4 + m25*ds5
          dq3 = m31*ds1 + m32*ds2 + m33*ds3 + m34*ds4 + m35*ds5
          dq4 = m41*ds1 + m42*ds2 + m43*ds3 + m44*ds4 + m45*ds5
          dq5 = m51*ds1 + m52*ds2 + m53*ds3 + m54*ds4 + m55*ds5
          if (neqns>=6) dq6 = ds6
          if (neqns>=7) dq7 = ds7
          if (neqns>=8) dq8 = ds8
          if (neqns>=9) dq9 = ds9
          if (neqns>=10)dq10= ds10
          if (neqns>=11)dq11= ds11
          if (neqns>=12)dq12= ds12

          FLUX_ALN(ii,jj,1) = 0.5d0*(f1l+f1r)
          FLUX_ALN(ii,jj,2) = 0.5d0*(f2l+f2r)
          FLUX_ALN(ii,jj,3) = 0.5d0*(f3l+f3r)
          FLUX_ALN(ii,jj,4) = 0.5d0*(f4l+f4r)
          FLUX_ALN(ii,jj,5) = 0.5d0*(f5l+f5r)
          if (neqns>=6) FLUX_ALN(ii,jj,6) = 0.5d0*(f6l+f6r)
          if (neqns>=7) FLUX_ALN(ii,jj,7) = 0.5d0*(f7l+f7r)
          if (neqns>=8) FLUX_ALN(ii,jj,8) = 0.5d0*(f8l+f8r)
          if (neqns>=9) FLUX_ALN(ii,jj,9) = 0.5d0*(f9l+f9r)
          if (neqns>=10)FLUX_ALN(ii,jj,10)= 0.5d0*(f10l+f10r)
          if (neqns>=11)FLUX_ALN(ii,jj,11)= 0.5d0*(f11l+f11r)
          if (neqns>=12)FLUX_ALN(ii,jj,12)= 0.5d0*(f12l+f12r)


          DLUX_ALN(ii,jj,1) =-0.5d0*dq1
          DLUX_ALN(ii,jj,2) =-0.5d0*dq2
          DLUX_ALN(ii,jj,3) =-0.5d0*dq3
          DLUX_ALN(ii,jj,4) =-0.5d0*dq4
          DLUX_ALN(ii,jj,5) =-0.5d0*dq5
          if (neqns>=6) DLUX_ALN(ii,jj,6) =-0.5d0*dq6
          if (neqns>=7) DLUX_ALN(ii,jj,7) =-0.5d0*dq7
          if (neqns>=8) DLUX_ALN(ii,jj,8) =-0.5d0*dq8
          if (neqns>=9) DLUX_ALN(ii,jj,9) =-0.5d0*dq9
          if (neqns>=10)DLUX_ALN(ii,jj,10)=-0.5d0*dq10
          if (neqns>=11)DLUX_ALN(ii,jj,11)=-0.5d0*dq11
          if (neqns>=12)DLUX_ALN(ii,jj,12)=-0.5d0*dq12
        enddo
      enddo

      do jj=1,jjmax
        do ii=1,iimax
          rhb1 =-(DLUX_ALN(ii,jj+1,1)-DLUX_ALN(ii,jj,1))
          rhb2 =-(DLUX_ALN(ii,jj+1,2)-DLUX_ALN(ii,jj,2))
          rhb3 =-(DLUX_ALN(ii,jj+1,3)-DLUX_ALN(ii,jj,3))
          rhb4 =-(DLUX_ALN(ii,jj+1,4)-DLUX_ALN(ii,jj,4))
          rhb5 =-(DLUX_ALN(ii,jj+1,5)-DLUX_ALN(ii,jj,5))
          if (neqns>=6) rhb6 =-(DLUX_ALN(ii,jj+1, 6)-DLUX_ALN(ii,jj, 6))
          if (neqns>=7) rhb7 =-(DLUX_ALN(ii,jj+1, 7)-DLUX_ALN(ii,jj, 7))
          if (neqns>=8) rhb8 =-(DLUX_ALN(ii,jj+1, 8)-DLUX_ALN(ii,jj, 8))
          if (neqns>=9) rhb9 =-(DLUX_ALN(ii,jj+1, 9)-DLUX_ALN(ii,jj, 9))
          if (neqns>=10)rhb10=-(DLUX_ALN(ii,jj+1,10)-DLUX_ALN(ii,jj,10))
          if (neqns>=11)rhb11=-(DLUX_ALN(ii,jj+1,11)-DLUX_ALN(ii,jj,11))
          if (neqns>=12)rhb12=-(DLUX_ALN(ii,jj+1,12)-DLUX_ALN(ii,jj,12))

          u = QQ_ALN(ii,jj,2)
          v = QQ_ALN(ii,jj,3)
          w = QQ_ALN(ii,jj,4)
          if (neqns>=6) ys1 = QQ_ALN(ii,jj,6)
          if (neqns>=7) ys2 = QQ_ALN(ii,jj,7)
          if (neqns>=8) ys3 = QQ_ALN(ii,jj,8)
          if (neqns>=9) ys4 = QQ_ALN(ii,jj,9)
          if (neqns>=10)ys5 = QQ_ALN(ii,jj,10)
          if (neqns>=11)ys6 = QQ_ALN(ii,jj,11)
          if (neqns>=12)ys7 = QQ_ALN(ii,jj,12)

          rh = DNST_ALN(ii,jj)
          hs = HTPS_ALN(ii,jj)

          rht = DRTM_ALN(ii,jj)
          rhp = DRPR_ALN(ii,jj)
          hst = SHCP_ALN(ii,jj)*CPREF
          hsp = DHPR_ALN(ii,jj)

          c2 = rh*hst/(rh*rhp*hst+rht*(1.0d0-rh*hsp))
          c  = dsqrt(c2)

          vp = dmin1(VPRE_ALN(ii,jj),c)
          vp2= vp*vp

          qq  = u*u+v*v+w*w

          ht = hs+0.5d0*qq

          alph = 1.0d0-rh*hsp

          rhc1 = rhb1+alph*(c2/vp2-1.0d0)/(rh*hst)*rhb5
          rhc2 = rhb2
          rhc3 = rhb3
          rhc4 = rhb4
          rhc5 = c2/vp2*rhb5
          if (neqns>=6) rhc6 = rhb6
          if (neqns>=7) rhc7 = rhb7
          if (neqns>=8) rhc8 = rhb8
          if (neqns>=9) rhc9 = rhb9
          if (neqns>=10)rhc10= rhb10
          if (neqns>=11)rhc11= rhb11
          if (neqns>=12)rhc12= rhb12

          rhb1 =   rht*rhc1+rhp*rhc5
          rhb2 = u*rht*rhc1+rh*rhc2+u*rhp*rhc5
          rhb3 = v*rht*rhc1+rh*rhc3+v*rhp*rhc5
          rhb4 = w*rht*rhc1+rh*rhc4+w*rhp*rhc5
          rhb5 = (rht*ht+rh*hst)*rhc1+u*rh*rhc2+v*rh*rhc3+w*rh*rhc4+(ht*rhp-alph)*rhc5
          if (neqns>=6) rhb6 = ys1*rht*rhc1+ys1*rhp*rhc5+rh*rhc6
          if (neqns>=7) rhb7 = ys2*rht*rhc1+ys2*rhp*rhc5+rh*rhc7
          if (neqns>=8) rhb8 = ys3*rht*rhc1+ys3*rhp*rhc5+rh*rhc8
          if (neqns>=9) rhb9 = ys4*rht*rhc1+ys4*rhp*rhc5+rh*rhc9
          if (neqns>=10)rhb10= ys5*rht*rhc1+ys5*rhp*rhc5+rh*rhc10
          if (neqns>=11)rhb11= ys6*rht*rhc1+ys6*rhp*rhc5+rh*rhc11
          if (neqns>=12)rhb12= ys7*rht*rhc1+ys7*rhp*rhc5+rh*rhc12

          DFLX_ALN(ii,jj,1) =-(FLUX_ALN(ii,jj+1,1)-FLUX_ALN(ii,jj,1)) + rhb1
          DFLX_ALN(ii,jj,2) =-(FLUX_ALN(ii,jj+1,2)-FLUX_ALN(ii,jj,2)) + rhb2
          DFLX_ALN(ii,jj,3) =-(FLUX_ALN(ii,jj+1,3)-FLUX_ALN(ii,jj,3)) + rhb3
          DFLX_ALN(ii,jj,4) =-(FLUX_ALN(ii,jj+1,4)-FLUX_ALN(ii,jj,4)) + rhb4
          DFLX_ALN(ii,jj,5) =-(FLUX_ALN(ii,jj+1,5)-FLUX_ALN(ii,jj,5)) + rhb5
          if (neqns>=6) DFLX_ALN(ii,jj,6) =-(FLUX_ALN(ii,jj+1, 6)-FLUX_ALN(ii,jj, 6)) + rhb6
          if (neqns>=7) DFLX_ALN(ii,jj,7) =-(FLUX_ALN(ii,jj+1, 7)-FLUX_ALN(ii,jj, 7)) + rhb7
          if (neqns>=8) DFLX_ALN(ii,jj,8) =-(FLUX_ALN(ii,jj+1, 8)-FLUX_ALN(ii,jj, 8)) + rhb8
          if (neqns>=9) DFLX_ALN(ii,jj,9) =-(FLUX_ALN(ii,jj+1, 9)-FLUX_ALN(ii,jj, 9)) + rhb9
          if (neqns>=10)DFLX_ALN(ii,jj,10)=-(FLUX_ALN(ii,jj+1,10)-FLUX_ALN(ii,jj,10)) + rhb10
          if (neqns>=11)DFLX_ALN(ii,jj,11)=-(FLUX_ALN(ii,jj+1,11)-FLUX_ALN(ii,jj,11)) + rhb11
          if (neqns>=12)DFLX_ALN(ii,jj,12)=-(FLUX_ALN(ii,jj+1,12)-FLUX_ALN(ii,jj,12)) + rhb12
        enddo
      enddo

      deallocate(DTMP)
      deallocate(THTA)

    end subroutine calc_fvs

  end subroutine flux_riemann_fvs_real_gas_extension

end module mdl_flux_riemann_fvs_real_gas_extension
