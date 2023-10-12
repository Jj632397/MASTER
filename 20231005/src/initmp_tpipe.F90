module mdl_initmp
  use mdl_param
  use mdl_refs
  use mdl_look_up_table
  use mdl_table
  implicit none

  !main passage
  real(8), public, save :: RHUNI_MP
  real(8), public, save :: UUUNI_MP
  real(8), public, save :: PPUNI_MP
  real(8), public, save :: TTUNI_MP
  real(8), public, save :: VSUNI_MP
  real(8), public, save :: CPUNI_MP
  real(8), public, save :: CCUNI_MP
  real(8), public, save :: MACHN_MP

  !sub passage
  real(8), public, save :: RHUNI_SP
  real(8), public, save :: UUUNI_SP
  real(8), public, save :: PPUNI_SP
  real(8), public, save :: TTUNI_SP
  real(8), public, save :: VSUNI_SP
  real(8), public, save :: CPUNI_SP
  real(8), public, save :: CCUNI_SP
  real(8), public, save :: MACHN_SP

  !passive scalar flowing from sub passage
  real(8), public, save :: SCALR_SP

  !back pressure
  real(8), public, save :: PPBCK

contains

  subroutine initmp
    integer :: ivalid
    real(8) :: rh,rht,rhp,hst,hsp,c2
    real(8) :: Q_MP, Q_SP, WDTH_MP,S_MP,UAVE_MP,UAVE_SP,RHTEMP_MP,RHTEMP_SP,VSTEMP_MP  !FURUSAWA

    SCALR_SP = 1.0d0
    
    PPBCK = 30.0d6   !back pressure

    TTUNI_MP = 703.12D0  ! 400 condition. inflow temperature of main passage
!    TTUNI_MP = 677.61D0  ! 380 condition. inflow temperature of main passage
    TTUNI_SP =  20.0D0+273.15D0 !inflow temperature of sub passage

    Q_MP =  22.0D0 *1.0D-3/60 ! [g/min] -> [kg/s]
    Q_SP =  Q_MP/4.0D0 

    call look_up_table_calc_func_single(TBL_TP_to_DNST,TTUNI_MP,1.0d0,PPBCK,1.0d0,RHTEMP_MP,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_DNST,TTUNI_SP,1.0d0,PPBCK,1.0d0,RHTEMP_SP,1.0d0,ivalid)

    WDTH_MP = 1.3D0 *1.0D-3 ! [mm] -> [m]
    S_MP  = WDTH_MP/2.0D0*WDTH_MP/2.0D0*4.0D0*DATAN(1.0D0) ! AREA for Pipe
    UAVE_MP = Q_MP/(S_MP*RHTEMP_MP)
    UAVE_SP = Q_SP/(S_MP*RHTEMP_SP)

    UUUNI_MP =  UAVE_MP *1.5D0 !inflow velocity of main passage
    UUUNI_SP =  UAVE_SP *1.5D0 !inflow velocity of sub passage

!    UUUNI_MP =  30.0d0  !inflow velocity of main passage
!    UUUNI_SP =   1.0d0  !inflow velocity of sub passage

    call look_up_table_calc_func_single(TBL_TP_to_VSSV,TTUNI_MP,1.0d0,PPBCK,1.0d0,VSTEMP_MP,1.0d0,ivalid)

    RYNLD = RHTEMP_MP*UUUNI_MP*WDTH_MP/VSTEMP_MP        !Reynold's number of main passage

!    RYNLD = 1.0d4       !Reynold's number of main passage

    !main passage
    call look_up_table_calc_func_single(TBL_TP_to_DNST,TTUNI_MP,1.0d0,PPBCK,1.0d0,RHUNI_MP,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_SHCP,TTUNI_MP,1.0d0,PPBCK,1.0d0,CPUNI_MP,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_VSSV,TTUNI_MP,1.0d0,PPBCK,1.0d0,VSUNI_MP,1.0d0,ivalid)

    call look_up_table_calc_func_single(TBL_TP_to_DRTM,TTUNI_MP,1.0d0,PPBCK,1.0d0,rht,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_DRPR,TTUNI_MP,1.0d0,PPBCK,1.0d0,rhp,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_DHPR,TTUNI_MP,1.0d0,PPBCK,1.0d0,hsp,1.0d0,ivalid)

    rh = RHUNI_MP
    hst= CPUNI_MP

    c2 = rh*hst/(rh*rhp*hst+rht*(1.0d0-rh*hsp))
    CCUNI_MP = dsqrt(c2)

    MACHN_MP = UUUNI_MP/CCUNI_MP


    !sub passage
    call look_up_table_calc_func_single(TBL_TP_to_DNST,TTUNI_SP,1.0d0,PPBCK,1.0d0,RHUNI_SP,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_SHCP,TTUNI_SP,1.0d0,PPBCK,1.0d0,CPUNI_SP,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_VSSV,TTUNI_SP,1.0d0,PPBCK,1.0d0,VSUNI_SP,1.0d0,ivalid)

    call look_up_table_calc_func_single(TBL_TP_to_DRTM,TTUNI_SP,1.0d0,PPBCK,1.0d0,rht,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_DRPR,TTUNI_SP,1.0d0,PPBCK,1.0d0,rhp,1.0d0,ivalid)
    call look_up_table_calc_func_single(TBL_TP_to_DHPR,TTUNI_SP,1.0d0,PPBCK,1.0d0,hsp,1.0d0,ivalid)

    rh = RHUNI_SP
    hst= CPUNI_SP

    c2 = rh*hst/(rh*rhp*hst+rht*(1.0d0-rh*hsp))
    CCUNI_SP = dsqrt(c2)

    MACHN_SP = UUUNI_SP/CCUNI_SP

    !Velocity(Time) scale for precodition
    MACHP = MACHN_MP
    UUPRE = CCUNI_MP*MACHP

    MACHN = MACHN_MP
    !------------
    RHREF = RHUNI_MP
    UUREF = UUUNI_MP
    CPREF = CPUNI_MP
    VSREF = VSUNI_MP
    !-----------
    TTREF = UUREF*UUREF
    PPREF = RHREF*UUREF*UUREF
    XYREF = RYNLD*VSREF/(RHREF*UUREF)
    CKREF = VSREF*CPREF/prntl
    DTREF = XYREF/UUREF
    DFREF = UUREF*XYREF
  end subroutine initmp

end module mdl_initmp
