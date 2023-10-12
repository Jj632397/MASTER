module mdl_block
  use mdl_param
  implicit none

  type block
    !Exclude halo region
    integer :: imax,jmax,kmax
    integer :: ista,jsta,ksta
    integer :: iend,jend,kend

    !Include halo region
    integer :: imax_wh,jmax_wh,kmax_wh
    integer :: ista_wh,jsta_wh,ksta_wh
    integer :: iend_wh,jend_wh,kend_wh

    integer, dimension(3) :: iperiodic

    !Node values----------------------------------------------------------------
    real(8), pointer, dimension(:,:,:)   :: XCRD_NODE,YCRD_NODE,ZCRD_NODE !x y z

    !Cell center values---------------------------------------------------------
    !Conserved quantities
    real(8), pointer, dimension(:,:,:,:) :: QCS1 !n   step
    real(8), pointer, dimension(:,:,:,:) :: QCS2 !n-1 step
    real(8), pointer, dimension(:,:,:,:) :: QCS3 !n-2 step

    !Residuals
    real(8), pointer, dimension(:,:,:,:) :: RHSS

    !Arbitrary quantities for output (neqns*3)
    real(8), pointer, dimension(:,:,:,:) :: QOUT

    real(8), pointer, dimension(:,:,:)   :: XCRD,YCRD,ZCRD !x y z
    real(8), pointer, dimension(:,:,:)   :: DTLC !dt
    real(8), pointer, dimension(:,:,:)   :: ULCT,VLCT,WLCT !u v w
    real(8), pointer, dimension(:,:,:)   :: DNST !density
    real(8), pointer, dimension(:,:,:)   :: PRSS !pressure
    real(8), pointer, dimension(:,:,:)   :: TMPR !temperature
    real(8), pointer, dimension(:,:,:)   :: HTPS !specific static enthalpy
    real(8), pointer, dimension(:,:,:)   :: VSSV !shear viscous coefficient
    real(8), pointer, dimension(:,:,:)   :: VSBV !bulk viscous coefficient
    real(8), pointer, dimension(:,:,:)   :: CKPV !heat transfer coefficient
    real(8), pointer, dimension(:,:,:)   :: VSST !turbulence shear viscous coefficient
    real(8), pointer, dimension(:,:,:)   :: VSBT !turbulence bulk viscous coefficient
    real(8), pointer, dimension(:,:,:)   :: CKPT !turbulence heat transfer coefficient
    real(8), pointer, dimension(:,:,:)   :: SHCP !Specific heat at constant pressure
    real(8), pointer, dimension(:,:,:)   :: CSOS !speed of sound
    real(8), pointer, dimension(:,:,:)   :: DIEL !dielectric constant
    real(8), pointer, dimension(:,:,:)   :: KNET !reaction rate

    real(8), pointer, dimension(:,:,:)   :: DHPR !h,p
    real(8), pointer, dimension(:,:,:)   :: DRTM !rh,T
    real(8), pointer, dimension(:,:,:)   :: DRPR !rh,p

    real(8), pointer, dimension(:,:,:)   :: RJCB

    real(8), pointer, dimension(:,:,:)   :: GRSC

    real(8), pointer, dimension(:,:,:)   :: VPRE !Velocity(Time) scale for preconditioning

    real(8), pointer, dimension(:,:,:,:) :: YSPC !scalar
    real(8), pointer, dimension(:,:,:,:) :: DYSX
    real(8), pointer, dimension(:,:,:,:) :: DYSY
    real(8), pointer, dimension(:,:,:,:) :: DYSZ

    !-----------------------------------------
    real(8), pointer, dimension(:,:,:) :: YNKK !distance_from_wall
    real(8), pointer, dimension(:,:,:) :: HFWL !wall dump function
    !-----------------------------------------

    real(8), pointer, dimension(:,:,:)   :: DXSX,DXSY,DXSZ !xsi,x xsi,y xsi,z
    real(8), pointer, dimension(:,:,:)   :: DETX,DETY,DETZ !eta,x eta,y eta,z
    real(8), pointer, dimension(:,:,:)   :: DZTX,DZTY,DZTZ !zta,x zta,y zta,z

    real(8), pointer, dimension(:,:,:)   :: DDNX,DDNY,DDNZ !rh,x rh,y rh,z
    real(8), pointer, dimension(:,:,:)   :: DULX,DULY,DULZ !u,x u,y u,z
    real(8), pointer, dimension(:,:,:)   :: DVLX,DVLY,DVLZ !v,x v,y v,z
    real(8), pointer, dimension(:,:,:)   :: DWLX,DWLY,DWLZ !w,x w,y w,z
    real(8), pointer, dimension(:,:,:)   :: DPRX,DPRY,DPRZ !w,x w,y w,z
    real(8), pointer, dimension(:,:,:)   :: DTMX,DTMY,DTMZ !T,x T,y T,z

    !Cell face values-----------------------------------------------------------
    real(8), pointer, dimension(:,:,:)   :: DXSX_FACE,DXSY_FACE,DXSZ_FACE !xsi,x xsi,y xsi,z
    real(8), pointer, dimension(:,:,:)   :: DETX_FACE,DETY_FACE,DETZ_FACE !eta,x eta,y eta,z
    real(8), pointer, dimension(:,:,:)   :: DZTX_FACE,DZTY_FACE,DZTZ_FACE !zta,x zta,y zta,z

    real(8), pointer, dimension(:,:,:)   :: XCRD_FACE_ETAZTA,YCRD_FACE_ETAZTA,ZCRD_FACE_ETAZTA !x y z (face zta-eta)
    real(8), pointer, dimension(:,:,:)   :: XCRD_FACE_ZTAXSI,YCRD_FACE_ZTAXSI,ZCRD_FACE_ZTAXSI !x y z (face zta-xsi)
    real(8), pointer, dimension(:,:,:)   :: XCRD_FACE_XSIETA,YCRD_FACE_XSIETA,ZCRD_FACE_XSIETA !x y z (face xsi-eta)

  end type block

contains

  subroutine block_allocate(BLK,ni,nj,nk)
    type(block), intent(inout) :: BLK
    integer, intent(in) :: ni,nj,nk

    integer :: nl
    integer :: ni_wh,nj_wh,nk_wh

    nl = neqns

    !-----------------------

    BLK%imax = ni
    BLK%jmax = nj
    BLK%kmax = nk

    BLK%ista = halo_width+1
    BLK%iend = halo_width+ni

    BLK%jsta = halo_width+1
    BLK%jend = halo_width+nj

    BLK%ksta = halo_width+1
    BLK%kend = halo_width+nk
    if (i2dimension>0) then
      BLK%ksta = 1
      BLK%kend = 1
    endif

    !-----------------------

    BLK%imax_wh = ni+halo_width*2
    BLK%jmax_wh = nj+halo_width*2
    BLK%kmax_wh = nk+halo_width*2
    if (i2dimension>0) then
      BLK%kmax_wh = 1
    endif

    BLK%ista_wh = 1
    BLK%iend_wh = ni+halo_width*2

    BLK%jsta_wh = 1
    BLK%jend_wh = nj+halo_width*2

    BLK%ksta_wh = 1
    BLK%kend_wh = nk+halo_width*2
    if (i2dimension>0) then
      BLK%ksta_wh = 1
      BLK%kend_wh = 1
    endif

    BLK%iperiodic(1) = 0
    BLK%iperiodic(2) = 0
    BLK%iperiodic(3) = 0

    ni_wh = ni+halo_width*2
    nj_wh = nj+halo_width*2
    nk_wh = nk+halo_width*2
    if (i2dimension>0) then
      nk_wh = 1
    endif

    !Node values----------------------------------------------------------------
    allocate(BLK%XCRD_NODE(ni_wh+1,nj_wh+1,nk_wh+1))
    allocate(BLK%YCRD_NODE(ni_wh+1,nj_wh+1,nk_wh+1))
    allocate(BLK%ZCRD_NODE(ni_wh+1,nj_wh+1,nk_wh+1))

    !Cell center values---------------------------------------------------------
    allocate(BLK%QCS1(ni_wh,nj_wh,nk_wh,nl))
    allocate(BLK%QCS2(ni_wh,nj_wh,nk_wh,nl))
    allocate(BLK%QCS3(ni_wh,nj_wh,nk_wh,nl))
    !---------------------------------------
    allocate(BLK%RHSS(ni_wh,nj_wh,nk_wh,nl))
    !---------------------------------------
    allocate(BLK%QOUT(ni_wh,nj_wh,nk_wh,nl*3))
    !------------------------------------
    allocate(BLK%DTLC(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%XCRD(ni_wh,nj_wh,nk_wh))
    allocate(BLK%YCRD(ni_wh,nj_wh,nk_wh))
    allocate(BLK%ZCRD(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%ULCT(ni_wh,nj_wh,nk_wh))
    allocate(BLK%VLCT(ni_wh,nj_wh,nk_wh))
    allocate(BLK%WLCT(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%DNST(ni_wh,nj_wh,nk_wh))
    allocate(BLK%PRSS(ni_wh,nj_wh,nk_wh))
    allocate(BLK%TMPR(ni_wh,nj_wh,nk_wh))
    allocate(BLK%HTPS(ni_wh,nj_wh,nk_wh))
    allocate(BLK%VSSV(ni_wh,nj_wh,nk_wh))
    allocate(BLK%VSBV(ni_wh,nj_wh,nk_wh))
    allocate(BLK%CKPV(ni_wh,nj_wh,nk_wh))
    allocate(BLK%VSST(ni_wh,nj_wh,nk_wh))
    allocate(BLK%VSBT(ni_wh,nj_wh,nk_wh))
    allocate(BLK%CKPT(ni_wh,nj_wh,nk_wh))
    allocate(BLK%SHCP(ni_wh,nj_wh,nk_wh))
    allocate(BLK%CSOS(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DIEL(ni_wh,nj_wh,nk_wh))
    allocate(BLK%KNET(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%DHPR(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DRTM(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DRPR(ni_wh,nj_wh,nk_wh))

    allocate(BLK%RJCB(ni_wh,nj_wh,nk_wh))

    allocate(BLK%GRSC(ni_wh,nj_wh,nk_wh))

    allocate(BLK%VPRE(ni_wh,nj_wh,nk_wh))

    if (nl-5>0) then
      allocate(BLK%YSPC(ni_wh,nj_wh,nk_wh,nl-5))
      allocate(BLK%DYSX(ni_wh,nj_wh,nk_wh,nl-5))
      allocate(BLK%DYSY(ni_wh,nj_wh,nk_wh,nl-5))
      allocate(BLK%DYSZ(ni_wh,nj_wh,nk_wh,nl-5))
    endif

    allocate(BLK%YNKK(ni_wh,nj_wh,nk_wh))
    allocate(BLK%HFWL(ni_wh,nj_wh,nk_wh))

    allocate(BLK%DXSX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DXSY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DXSZ(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DETX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DETY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DETZ(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DZTX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DZTY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DZTZ(ni_wh,nj_wh,nk_wh))

    allocate(BLK%DDNX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DDNY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DDNZ(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%DULX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DULY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DULZ(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%DVLX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DVLY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DVLZ(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%DWLX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DWLY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DWLZ(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%DPRX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DPRY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DPRZ(ni_wh,nj_wh,nk_wh))
    !------------------------------------
    allocate(BLK%DTMX(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DTMY(ni_wh,nj_wh,nk_wh))
    allocate(BLK%DTMZ(ni_wh,nj_wh,nk_wh))
    !------------------------------------

    !Cell face values-----------------------------------------------------------
    allocate(BLK%DXSX_FACE(ni_wh+1,nj_wh,nk_wh))
    allocate(BLK%DXSY_FACE(ni_wh+1,nj_wh,nk_wh))
    allocate(BLK%DXSZ_FACE(ni_wh+1,nj_wh,nk_wh))
    !-------------------------------------------
    allocate(BLK%DETX_FACE(ni_wh,nj_wh+1,nk_wh))
    allocate(BLK%DETY_FACE(ni_wh,nj_wh+1,nk_wh))
    allocate(BLK%DETZ_FACE(ni_wh,nj_wh+1,nk_wh))
    !-------------------------------------------
    allocate(BLK%DZTX_FACE(ni_wh,nj_wh,nk_wh+1))
    allocate(BLK%DZTY_FACE(ni_wh,nj_wh,nk_wh+1))
    allocate(BLK%DZTZ_FACE(ni_wh,nj_wh,nk_wh+1))

    allocate(BLK%XCRD_FACE_ETAZTA(ni_wh+1,nj_wh,nk_wh))
    allocate(BLK%YCRD_FACE_ETAZTA(ni_wh+1,nj_wh,nk_wh))
    allocate(BLK%ZCRD_FACE_ETAZTA(ni_wh+1,nj_wh,nk_wh))
    !--------------------------------------------------
    allocate(BLK%XCRD_FACE_ZTAXSI(ni_wh,nj_wh+1,nk_wh))
    allocate(BLK%YCRD_FACE_ZTAXSI(ni_wh,nj_wh+1,nk_wh))
    allocate(BLK%ZCRD_FACE_ZTAXSI(ni_wh,nj_wh+1,nk_wh))
    !--------------------------------------------------
    allocate(BLK%XCRD_FACE_XSIETA(ni_wh,nj_wh,nk_wh+1))
    allocate(BLK%YCRD_FACE_XSIETA(ni_wh,nj_wh,nk_wh+1))
    allocate(BLK%ZCRD_FACE_XSIETA(ni_wh,nj_wh,nk_wh+1))

    BLK%XCRD_NODE = 0.0d0
    BLK%YCRD_NODE = 0.0d0
    BLK%ZCRD_NODE = 0.0d0

    BLK%QCS1 = 0.0d0
    BLK%QCS2 = 0.0d0
    BLK%QCS3 = 0.0d0
    !---------------
    BLK%RHSS = 0.0d0
    !---------------
    BLK%QOUT = 0.0d0
    !---------------
    BLK%DTLC = 0.0d0
    !---------------
    BLK%XCRD = 0.0d0
    BLK%YCRD = 0.0d0
    BLK%ZCRD = 0.0d0
    !---------------
    BLK%ULCT = 0.0d0
    BLK%VLCT = 0.0d0
    BLK%WLCT = 0.0d0
    !---------------
    BLK%DNST = 0.0d0
    BLK%PRSS = 0.0d0
    BLK%TMPR = 0.0d0
    BLK%HTPS = 0.0d0
    BLK%VSSV = 0.0d0
    BLK%VSBV = 0.0d0
    BLK%CKPV = 0.0d0
    BLK%VSST = 0.0d0
    BLK%VSBT = 0.0d0
    BLK%CKPT = 0.0d0
    BLK%SHCP = 0.0d0
    BLK%CSOS = 0.0d0
    BLK%DIEL = 0.0d0
    BLK%KNET = 0.0d0
    !---------------
    BLK%DHPR = 0.0d0
    BLK%DRTM = 0.0d0
    BLK%DRPR = 0.0d0
    !---------------
    BLK%RJCB = 0.0d0
    BLK%GRSC = 0.0d0
    !---------------
    BLK%VPRE = 0.0d0
    !---------------
    if (nl-5>0) then
      BLK%YSPC = 0.0d0
      BLK%DYSX = 0.0d0
      BLK%DYSY = 0.0d0
      BLK%DYSZ = 0.0d0
    endif
    !---------------
    BLK%YNKK = 0.0d0
    BLK%HFWL = 0.0d0
    !---------------
    BLK%DXSX = 0.0d0
    BLK%DXSY = 0.0d0
    BLK%DXSZ = 0.0d0
    BLK%DETX = 0.0d0
    BLK%DETY = 0.0d0
    BLK%DETZ = 0.0d0
    BLK%DZTX = 0.0d0
    BLK%DZTY = 0.0d0
    BLK%DZTZ = 0.0d0
    !---------------
    BLK%DDNX = 0.0d0
    BLK%DDNY = 0.0d0
    BLK%DDNZ = 0.0d0
    !---------------
    BLK%DULX = 0.0d0
    BLK%DULY = 0.0d0
    BLK%DULZ = 0.0d0
    !---------------
    BLK%DVLX = 0.0d0
    BLK%DVLY = 0.0d0
    BLK%DVLZ = 0.0d0
    !---------------
    BLK%DWLX = 0.0d0
    BLK%DWLY = 0.0d0
    BLK%DWLZ = 0.0d0
    !---------------
    BLK%DPRX = 0.0d0
    BLK%DPRY = 0.0d0
    BLK%DPRZ = 0.0d0
    !---------------
    BLK%DTMX = 0.0d0
    BLK%DTMY = 0.0d0
    BLK%DTMZ = 0.0d0
    !---------------
    BLK%DXSX_FACE = 0.0d0
    BLK%DXSY_FACE = 0.0d0
    BLK%DXSZ_FACE = 0.0d0
    !--------------------
    BLK%DETX_FACE = 0.0d0
    BLK%DETY_FACE = 0.0d0
    BLK%DETZ_FACE = 0.0d0
    !--------------------
    BLK%DZTX_FACE = 0.0d0
    BLK%DZTY_FACE = 0.0d0
    BLK%DZTZ_FACE = 0.0d0
    !--------------------
    BLK%XCRD_FACE_ETAZTA = 0.0d0
    BLK%YCRD_FACE_ETAZTA = 0.0d0
    BLK%ZCRD_FACE_ETAZTA = 0.0d0
    !---------------------------
    BLK%XCRD_FACE_ZTAXSI = 0.0d0
    BLK%YCRD_FACE_ZTAXSI = 0.0d0
    BLK%ZCRD_FACE_ZTAXSI = 0.0d0
    !---------------------------
    BLK%XCRD_FACE_XSIETA = 0.0d0
    BLK%YCRD_FACE_XSIETA = 0.0d0
    BLK%ZCRD_FACE_XSIETA = 0.0d0
  end subroutine block_allocate

  subroutine block_deallocate(BLK)
    type(block), intent(inout) :: BLK

    deallocate(BLK%XCRD_NODE)
    deallocate(BLK%YCRD_NODE)
    deallocate(BLK%ZCRD_NODE)

    deallocate(BLK%QCS1)
    deallocate(BLK%QCS2)
    deallocate(BLK%QCS3)
    !-------------------
    deallocate(BLK%RHSS)
    !-------------------
    deallocate(BLK%QOUT)
    !-------------------
    deallocate(BLK%DTLC)
    !-------------------
    deallocate(BLK%XCRD)
    deallocate(BLK%YCRD)
    deallocate(BLK%ZCRD)
    !-------------------
    deallocate(BLK%ULCT)
    deallocate(BLK%VLCT)
    deallocate(BLK%WLCT)
    !-------------------
    deallocate(BLK%DNST)
    deallocate(BLK%PRSS)
    deallocate(BLK%TMPR)
    deallocate(BLK%HTPS)
    deallocate(BLK%VSSV)
    deallocate(BLK%VSBV)
    deallocate(BLK%CKPV)
    deallocate(BLK%VSST)
    deallocate(BLK%VSBT)
    deallocate(BLK%CKPT)
    deallocate(BLK%SHCP)
    deallocate(BLK%CSOS)
    deallocate(BLK%DIEL)
    deallocate(BLK%KNET)
    !-------------------
    deallocate(BLK%DHPR)
    deallocate(BLK%DRTM)
    deallocate(BLK%DRPR)
    deallocate(BLK%RJCB)
    deallocate(BLK%GRSC)
    deallocate(BLK%VPRE)
    !-------------------
    if (neqns-5>0) then
      deallocate(BLK%YSPC)
      deallocate(BLK%DYSX)
      deallocate(BLK%DYSY)
      deallocate(BLK%DYSZ)
    endif
    !-------------------
    deallocate(BLK%YNKK)
    deallocate(BLK%HFWL)
    !-------------------
    deallocate(BLK%DXSX)
    deallocate(BLK%DXSY)
    deallocate(BLK%DXSZ)
    deallocate(BLK%DETX)
    deallocate(BLK%DETY)
    deallocate(BLK%DETZ)
    deallocate(BLK%DZTX)
    deallocate(BLK%DZTY)
    deallocate(BLK%DZTZ)
    !-------------------
    deallocate(BLK%DDNX)
    deallocate(BLK%DDNY)
    deallocate(BLK%DDNZ)
    !-------------------
    deallocate(BLK%DULX)
    deallocate(BLK%DULY)
    deallocate(BLK%DULZ)
    !-------------------
    deallocate(BLK%DVLX)
    deallocate(BLK%DVLY)
    deallocate(BLK%DVLZ)
    !-------------------
    deallocate(BLK%DWLX)
    deallocate(BLK%DWLY)
    deallocate(BLK%DWLZ)
    !-------------------
    deallocate(BLK%DPRX)
    deallocate(BLK%DPRY)
    deallocate(BLK%DPRZ)
    !-------------------
    deallocate(BLK%DTMX)
    deallocate(BLK%DTMY)
    deallocate(BLK%DTMZ)
    !-------------------
    deallocate(BLK%DXSX_FACE)
    deallocate(BLK%DXSY_FACE)
    deallocate(BLK%DXSZ_FACE)
    !------------------------
    deallocate(BLK%DETX_FACE)
    deallocate(BLK%DETY_FACE)
    deallocate(BLK%DETZ_FACE)
    !------------------------
    deallocate(BLK%DZTX_FACE)
    deallocate(BLK%DZTY_FACE)
    deallocate(BLK%DZTZ_FACE)
    !------------------------
    deallocate(BLK%XCRD_FACE_ETAZTA)
    deallocate(BLK%YCRD_FACE_ETAZTA)
    deallocate(BLK%ZCRD_FACE_ETAZTA)
    !-------------------------------
    deallocate(BLK%XCRD_FACE_ZTAXSI)
    deallocate(BLK%YCRD_FACE_ZTAXSI)
    deallocate(BLK%ZCRD_FACE_ZTAXSI)
    !-------------------------------
    deallocate(BLK%XCRD_FACE_XSIETA)
    deallocate(BLK%YCRD_FACE_XSIETA)
    deallocate(BLK%ZCRD_FACE_XSIETA)
  end subroutine block_deallocate

end module mdl_block
