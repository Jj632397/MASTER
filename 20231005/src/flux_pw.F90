module mdl_flux_pw
  use mdl_mpisub_sbsp
  use mdl_decompo
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_fdm1
  use mdl_halo
  implicit none

  integer, private, parameter :: icompact = 1
  integer, private, parameter, dimension(5) :: isc_xsi = (/2,1,1,1,2/)
  integer, private, parameter, dimension(5) :: isc_eta = (/2,1,1,1,2/)
  integer, private, parameter, dimension(5) :: isc_zta = (/2,1,1,1,2/)

contains

  subroutine flux_pw(BLK,RHS)
    type(block), intent(inout) :: BLK
    real(8), intent(inout), dimension(:,:,:,:) :: RHS

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: iimax,jjmax

    real(8), allocatable, dimension(:,:)   :: DNST_ALN
    real(8), allocatable, dimension(:,:)   :: ULCT_ALN
    real(8), allocatable, dimension(:,:)   :: VLCT_ALN
    real(8), allocatable, dimension(:,:)   :: WLCT_ALN
    real(8), allocatable, dimension(:,:)   :: PRSS_ALN
    real(8), allocatable, dimension(:,:)   :: HTPS_ALN
    real(8), allocatable, dimension(:,:,:) :: YSPC_ALN

    real(8), allocatable, dimension(:,:)   :: DXSX_ALN
    real(8), allocatable, dimension(:,:)   :: DXSY_ALN
    real(8), allocatable, dimension(:,:)   :: DXSZ_ALN

    real(8), allocatable, dimension(:,:)   :: XCRD_ALN
    real(8), allocatable, dimension(:,:)   :: YCRD_ALN
    real(8), allocatable, dimension(:,:)   :: ZCRD_ALN

    real(8), allocatable, dimension(:,:,:) :: FLUX_ALN

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

      if (IFLAG_PRIM(1)==1) nhalo1 = halo_width
      if (IFLAG_PRIM(2)==1) nhalo2 = halo_width

      iimax = jmax*kmax
      jjmax = imax+nhalo1+nhalo2

      ioft = BLK%ista-1
      joft = BLK%jsta-1
      koft = BLK%ksta-1

      if (IFLAG_PRIM(1)==1) ioft = 0

      allocate(DNST_ALN(iimax,jjmax))
      allocate(ULCT_ALN(iimax,jjmax))
      allocate(VLCT_ALN(iimax,jjmax))
      allocate(WLCT_ALN(iimax,jjmax))
      allocate(PRSS_ALN(iimax,jjmax))
      allocate(HTPS_ALN(iimax,jjmax))
      if (neqns>5) then
        allocate(YSPC_ALN(iimax,jjmax,neqns-5))
      endif

      allocate(DXSX_ALN(iimax,jjmax))
      allocate(DXSY_ALN(iimax,jjmax))
      allocate(DXSZ_ALN(iimax,jjmax))

      allocate(XCRD_ALN(iimax,jjmax))
      allocate(YCRD_ALN(iimax,jjmax))
      allocate(ZCRD_ALN(iimax,jjmax))

      allocate(FLUX_ALN(iimax,jjmax,neqns))

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + jj
          j = joft + mod(ii-1,jmax)+1
          k = koft + (ii-1)/jmax+1
          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          ULCT_ALN(ii,jj) = BLK%ULCT(i,j,k)
          VLCT_ALN(ii,jj) = BLK%VLCT(i,j,k)
          WLCT_ALN(ii,jj) = BLK%WLCT(i,j,k)
          PRSS_ALN(ii,jj) = BLK%PRSS(i,j,k)
          HTPS_ALN(ii,jj) = BLK%HTPS(i,j,k)
          !--------------------------------
          DXSX_ALN(ii,jj) = BLK%DXSX(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DXSY(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DXSZ(i,j,k)
          !--------------------------------
          XCRD_ALN(ii,jj) = BLK%XCRD(i,j,k)
          YCRD_ALN(ii,jj) = BLK%YCRD(i,j,k)
          ZCRD_ALN(ii,jj) = BLK%ZCRD(i,j,k)
        enddo
      enddo

      do l=1,neqns-5
        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + jj
            j = joft + mod(ii-1,jmax)+1
            k = koft + (ii-1)/jmax+1
            YSPC_ALN(ii,jj,l) = BLK%YSPC(i,j,k,l)
          enddo
        enddo
      enddo

      call calc_flux_pw

      do l=1,neqns
        call calc_fdm_1st_align(FLUX_ALN(:,:,l),FLUX_ALN(:,:,l),0,icompact,isc_xsi,iimax,jjmax)
      enddo

      do l=1,neqns
        do jj=1,imax
          do ii=1,jmax*kmax
            i = BLK%ista-1 + jj
            j = BLK%jsta-1 + mod(ii-1,jmax)+1
            k = BLK%ksta-1 + (ii-1)/jmax+1
            RHS(i,j,k,l) = RHS(i,j,k,l) - FLUX_ALN(ii,nhalo1+jj,l)
          enddo
        enddo
      enddo

      deallocate(DNST_ALN)
      deallocate(ULCT_ALN)
      deallocate(VLCT_ALN)
      deallocate(WLCT_ALN)
      deallocate(PRSS_ALN)
      deallocate(HTPS_ALN)
      if (neqns>5) then
        deallocate(YSPC_ALN)
      endif
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      deallocate(XCRD_ALN)
      deallocate(YCRD_ALN)
      deallocate(ZCRD_ALN)
      deallocate(FLUX_ALN)
    end subroutine calc_xsi
!-------------------------------------------------------------------------------
    subroutine calc_eta
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft

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

      allocate(DNST_ALN(iimax,jjmax))
      allocate(ULCT_ALN(iimax,jjmax))
      allocate(VLCT_ALN(iimax,jjmax))
      allocate(WLCT_ALN(iimax,jjmax))
      allocate(PRSS_ALN(iimax,jjmax))
      allocate(HTPS_ALN(iimax,jjmax))
      if (neqns>5) then
        allocate(YSPC_ALN(iimax,jjmax,neqns-5))
      endif

      allocate(DXSX_ALN(iimax,jjmax))
      allocate(DXSY_ALN(iimax,jjmax))
      allocate(DXSZ_ALN(iimax,jjmax))

      allocate(XCRD_ALN(iimax,jjmax))
      allocate(YCRD_ALN(iimax,jjmax))
      allocate(ZCRD_ALN(iimax,jjmax))

      allocate(FLUX_ALN(iimax,jjmax,neqns))

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + jj
          k = koft + (ii-1)/imax+1
          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          ULCT_ALN(ii,jj) = BLK%ULCT(i,j,k)
          VLCT_ALN(ii,jj) = BLK%VLCT(i,j,k)
          WLCT_ALN(ii,jj) = BLK%WLCT(i,j,k)
          PRSS_ALN(ii,jj) = BLK%PRSS(i,j,k)
          HTPS_ALN(ii,jj) = BLK%HTPS(i,j,k)
          !--------------------------------
          DXSX_ALN(ii,jj) = BLK%DETX(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DETY(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DETZ(i,j,k)
          !--------------------------------
          XCRD_ALN(ii,jj) = BLK%XCRD(i,j,k)
          YCRD_ALN(ii,jj) = BLK%YCRD(i,j,k)
          ZCRD_ALN(ii,jj) = BLK%ZCRD(i,j,k)
        enddo
      enddo

      do l=1,neqns-5
        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + jj
            k = koft + (ii-1)/imax+1
            YSPC_ALN(ii,jj,l) = BLK%YSPC(i,j,k,l)
          enddo
        enddo
      enddo

      call calc_flux_pw

      do l=1,neqns
        call calc_fdm_1st_align(FLUX_ALN(:,:,l),FLUX_ALN(:,:,l),0,icompact,isc_eta,iimax,jjmax)
      enddo

      do l=1,neqns
        do jj=1,jmax
          do ii=1,imax*kmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + jj
            k = BLK%ksta-1 + (ii-1)/imax+1
            RHS(i,j,k,l) = RHS(i,j,k,l) - FLUX_ALN(ii,nhalo1+jj,l)
          enddo
        enddo
      enddo

      deallocate(DNST_ALN)
      deallocate(ULCT_ALN)
      deallocate(VLCT_ALN)
      deallocate(WLCT_ALN)
      deallocate(PRSS_ALN)
      deallocate(HTPS_ALN)
      if (neqns>5) then
        deallocate(YSPC_ALN)
      endif
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      deallocate(XCRD_ALN)
      deallocate(YCRD_ALN)
      deallocate(ZCRD_ALN)
      deallocate(FLUX_ALN)
    end subroutine calc_eta
!-------------------------------------------------------------------------------
    subroutine calc_zta
      integer :: i,j,k,l,ii,jj
      integer :: nhalo1,nhalo2
      integer :: ioft,joft,koft

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

      allocate(DNST_ALN(iimax,jjmax))
      allocate(ULCT_ALN(iimax,jjmax))
      allocate(VLCT_ALN(iimax,jjmax))
      allocate(WLCT_ALN(iimax,jjmax))
      allocate(PRSS_ALN(iimax,jjmax))
      allocate(HTPS_ALN(iimax,jjmax))
      if (neqns>5) then
        allocate(YSPC_ALN(iimax,jjmax,neqns-5))
      endif

      allocate(DXSX_ALN(iimax,jjmax))
      allocate(DXSY_ALN(iimax,jjmax))
      allocate(DXSZ_ALN(iimax,jjmax))

      allocate(XCRD_ALN(iimax,jjmax))
      allocate(YCRD_ALN(iimax,jjmax))
      allocate(ZCRD_ALN(iimax,jjmax))

      allocate(FLUX_ALN(iimax,jjmax,neqns))

      do jj=1,jjmax
        do ii=1,iimax
          i = ioft + mod(ii-1,imax)+1
          j = joft + (ii-1)/imax+1
          k = koft + jj
          DNST_ALN(ii,jj) = BLK%DNST(i,j,k)
          ULCT_ALN(ii,jj) = BLK%ULCT(i,j,k)
          VLCT_ALN(ii,jj) = BLK%VLCT(i,j,k)
          WLCT_ALN(ii,jj) = BLK%WLCT(i,j,k)
          PRSS_ALN(ii,jj) = BLK%PRSS(i,j,k)
          HTPS_ALN(ii,jj) = BLK%HTPS(i,j,k)
          !--------------------------------
          DXSX_ALN(ii,jj) = BLK%DZTX(i,j,k)
          DXSY_ALN(ii,jj) = BLK%DZTY(i,j,k)
          DXSZ_ALN(ii,jj) = BLK%DZTZ(i,j,k)
          !--------------------------------
          XCRD_ALN(ii,jj) = BLK%XCRD(i,j,k)
          YCRD_ALN(ii,jj) = BLK%YCRD(i,j,k)
          ZCRD_ALN(ii,jj) = BLK%ZCRD(i,j,k)
        enddo
      enddo

      do l=1,neqns-5
        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + (ii-1)/imax+1
            k = koft + jj
            YSPC_ALN(ii,jj,l) = BLK%YSPC(i,j,k,l)
          enddo
        enddo
      enddo

      call calc_flux_pw

      do l=1,neqns
        call calc_fdm_1st_align(FLUX_ALN(:,:,l),FLUX_ALN(:,:,l),0,icompact,isc_zta,iimax,jjmax)
      enddo

      do l=1,neqns
        do jj=1,kmax
          do ii=1,imax*jmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + (ii-1)/imax+1
            k = BLK%ksta-1 + jj
            RHS(i,j,k,l) = RHS(i,j,k,l) - FLUX_ALN(ii,jj+nhalo1,l)
          enddo
        enddo
      enddo

      deallocate(DNST_ALN)
      deallocate(ULCT_ALN)
      deallocate(VLCT_ALN)
      deallocate(WLCT_ALN)
      deallocate(PRSS_ALN)
      deallocate(HTPS_ALN)
      if (neqns>5) then
        deallocate(YSPC_ALN)
      endif
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      deallocate(XCRD_ALN)
      deallocate(YCRD_ALN)
      deallocate(ZCRD_ALN)
      deallocate(FLUX_ALN)
    end subroutine calc_zta
!-------------------------------------------------------------------------------
    subroutine calc_flux_pw
      integer :: ii,jj,l
      real(8) :: rh,u,v,w,p,hs,ht,uc,ys
      real(8) :: xsix,xsiy,xsiz

      real(8) :: pi
      real(8) :: nx,ny,nz,ox,oy,oz,rx,ry,rz
      real(8) :: a,b,g,ac

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

      do jj=1,jjmax
        do ii=1,iimax
          rh = DNST_ALN(ii,jj)
          u  = ULCT_ALN(ii,jj)
          v  = VLCT_ALN(ii,jj)
          w  = WLCT_ALN(ii,jj)
          p  = PRSS_ALN(ii,jj)
          hs = HTPS_ALN(ii,jj)
          !-------------------
          xsix = DXSX_ALN(ii,jj)
          xsiy = DXSY_ALN(ii,jj)
          xsiz = DXSZ_ALN(ii,jj)

          rx = XCRD_ALN(ii,jj)-frame_origin(1)
          ry = YCRD_ALN(ii,jj)-frame_origin(2)
          rz = ZCRD_ALN(ii,jj)-frame_origin(3)

          ht = hs+0.5d0*(u*u+v*v+w*w)

          a = oy*rz-oz*ry
          b = oz*rx-ox*rz
          g = ox*ry-oy*rx

          uc = u*xsix+v*xsiy+w*xsiz
          ac = a*xsix+b*xsiy+g*xsiz

          FLUX_ALN(ii,jj,1) = rh*(uc-ac)
          FLUX_ALN(ii,jj,2) = rh*(uc-ac)*u+xsix*p
          FLUX_ALN(ii,jj,3) = rh*(uc-ac)*v+xsiy*p
          FLUX_ALN(ii,jj,4) = rh*(uc-ac)*w+xsiz*p
          FLUX_ALN(ii,jj,5) = rh*(uc-ac)*ht+p*ac
        enddo
      enddo

      do l=1,neqns-5
        do jj=1,jjmax
          do ii=1,iimax
            rh = DNST_ALN(ii,jj)
            u  = ULCT_ALN(ii,jj)
            v  = VLCT_ALN(ii,jj)
            w  = WLCT_ALN(ii,jj)
            ys = YSPC_ALN(ii,jj,l)
            !-------------------
            xsix = DXSX_ALN(ii,jj)
            xsiy = DXSY_ALN(ii,jj)
            xsiz = DXSZ_ALN(ii,jj)

            rx = XCRD_ALN(ii,jj)-frame_origin(1)
            ry = YCRD_ALN(ii,jj)-frame_origin(2)
            rz = ZCRD_ALN(ii,jj)-frame_origin(3)

            a = oy*rz-oz*ry
            b = oz*rx-ox*rz
            g = ox*ry-oy*rx

            uc = u*xsix+v*xsiy+w*xsiz
            ac = a*xsix+b*xsiy+g*xsiz

            FLUX_ALN(ii,jj,l+5) = rh*(uc-ac)*ys
          enddo
        enddo
      enddo

    end subroutine calc_flux_pw

  end subroutine flux_pw

end module mdl_flux_pw
