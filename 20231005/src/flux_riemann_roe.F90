module mdl_flux_riemann_roe
  use mdl_mpisub_sbsp
  use mdl_decompo
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_muscl
  use mdl_weno
  use mdl_halo
  implicit none

!-------------------------------------------------------------------------------
!Roe's approximate Riemann solver with MUSCL or 5th-order WENO.
!Only for ideal gas use. Using cp_i and cr_i in param.F90 for gas property.
!-------------------------------------------------------------------------------

  integer, private, parameter :: isc_reconstruction = 1 !<1>:weno5 <else>:MUSCL

  integer, private, parameter :: iorder_muscl = 3

contains

  subroutine flux_riemann_roe(BLK,RHS)
    type(block), intent(inout) :: BLK
    real(8), intent(inout), dimension(:,:,:,:) :: RHS

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: iimax,jjmax

    !Cell center (temporal data array)-----------------
    real(8), allocatable, dimension(:,:) :: QQ_ALN_TEMP
    !Cell face (temporal data array)-------------------
    real(8), allocatable, dimension(:,:) :: QL_ALN_TEMP
    real(8), allocatable, dimension(:,:) :: QR_ALN_TEMP

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

    !Cell center
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
      !flux
      allocate(FLUX_ALN(iimax,jjmax+1,neqns))


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

      allocate(QQ_ALN_TEMP(iimax,jjmax))
      allocate(QL_ALN_TEMP(iimax,jjmax+1))
      allocate(QR_ALN_TEMP(iimax,jjmax+1))

      do l=1,neqns

        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + jj
            j = joft + mod(ii-1,jmax)+1
            k = koft + (ii-1)/jmax+1
            if (l==1) then
              QQ_ALN_TEMP(ii,jj) = BLK%DNST(i,j,k)
            else if (l==2) then
              QQ_ALN_TEMP(ii,jj) = BLK%ULCT(i,j,k)
            else if (l==3) then
              QQ_ALN_TEMP(ii,jj) = BLK%VLCT(i,j,k)
            else if (l==4) then
              QQ_ALN_TEMP(ii,jj) = BLK%WLCT(i,j,k)
            else if (l==5) then
              QQ_ALN_TEMP(ii,jj) = BLK%PRSS(i,j,k)
            else
              QQ_ALN_TEMP(ii,jj) = BLK%YSPC(i,j,k,l-5)
            endif
          enddo
        enddo

        if (isc_reconstruction==1) then
          call weno5_align(QQ_ALN_TEMP,QL_ALN_TEMP,QR_ALN_TEMP,iimax,jjmax)
        else
          call muscl_align(QQ_ALN_TEMP,QL_ALN_TEMP,QR_ALN_TEMP,iorder_muscl,iimax,jjmax)
        endif

        do jj=1,jjmax+1
          do ii=1,iimax
            QL_ALN(ii,jj,l) = QL_ALN_TEMP(ii,jj)
            QR_ALN(ii,jj,l) = QR_ALN_TEMP(ii,jj)
          enddo
        enddo

      enddo

      deallocate(QQ_ALN_TEMP)
      deallocate(QL_ALN_TEMP)
      deallocate(QR_ALN_TEMP)

      call calc_roe

      do l=1,neqns
        do jj=1,imax
          do ii=1,jmax*kmax
            i = BLK%ista-1 + jj
            j = BLK%jsta-1 + mod(ii-1,jmax)+1
            k = BLK%ksta-1 + (ii-1)/jmax+1
            RHS(i,j,k,l) = RHS(i,j,k,l) - (FLUX_ALN(ii,nhalo1+jj+1,l)-FLUX_ALN(ii,nhalo1+jj,l))
          enddo
        enddo
      enddo

      deallocate(QL_ALN)
      deallocate(QR_ALN)
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


      !Cell face----------------------------------
      !reconstructing values
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
      !flux
      allocate(FLUX_ALN(iimax,jjmax+1,neqns))


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

      allocate(QQ_ALN_TEMP(iimax,jjmax))
      allocate(QL_ALN_TEMP(iimax,jjmax+1))
      allocate(QR_ALN_TEMP(iimax,jjmax+1))

      do l=1,neqns

        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + jj
            k = koft + (ii-1)/imax+1
            if (l==1) then
              QQ_ALN_TEMP(ii,jj) = BLK%DNST(i,j,k)
            else if (l==2) then
              QQ_ALN_TEMP(ii,jj) = BLK%ULCT(i,j,k)
            else if (l==3) then
              QQ_ALN_TEMP(ii,jj) = BLK%VLCT(i,j,k)
            else if (l==4) then
              QQ_ALN_TEMP(ii,jj) = BLK%WLCT(i,j,k)
            else if (l==5) then
              QQ_ALN_TEMP(ii,jj) = BLK%PRSS(i,j,k)
            else
              QQ_ALN_TEMP(ii,jj) = BLK%YSPC(i,j,k,l-5)
            endif
          enddo
        enddo

        if (isc_reconstruction==1) then
          call weno5_align(QQ_ALN_TEMP,QL_ALN_TEMP,QR_ALN_TEMP,iimax,jjmax)
        else
          call muscl_align(QQ_ALN_TEMP,QL_ALN_TEMP,QR_ALN_TEMP,iorder_muscl,iimax,jjmax)
        endif

        do jj=1,jjmax+1
          do ii=1,iimax
            QL_ALN(ii,jj,l) = QL_ALN_TEMP(ii,jj)
            QR_ALN(ii,jj,l) = QR_ALN_TEMP(ii,jj)
          enddo
        enddo

      enddo

      deallocate(QQ_ALN_TEMP)
      deallocate(QL_ALN_TEMP)
      deallocate(QR_ALN_TEMP)

      call calc_roe

      do l=1,neqns
        do jj=1,jmax
          do ii=1,imax*kmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + jj
            k = BLK%ksta-1 + (ii-1)/imax+1
            RHS(i,j,k,l) = RHS(i,j,k,l) - (FLUX_ALN(ii,nhalo1+jj+1,l)-FLUX_ALN(ii,nhalo1+jj,l))
          enddo
        enddo
      enddo

      deallocate(QL_ALN)
      deallocate(QR_ALN)
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


      !Cell face----------------------------------
      !reconstructing values
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
      !flux
      allocate(FLUX_ALN(iimax,jjmax+1,neqns))


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

      allocate(QQ_ALN_TEMP(iimax,jjmax))
      allocate(QL_ALN_TEMP(iimax,jjmax+1))
      allocate(QR_ALN_TEMP(iimax,jjmax+1))

      do l=1,neqns

        do jj=1,jjmax
          do ii=1,iimax
            i = ioft + mod(ii-1,imax)+1
            j = joft + (ii-1)/imax+1
            k = koft + jj
            if (l==1) then
              QQ_ALN_TEMP(ii,jj) = BLK%DNST(i,j,k)
            else if (l==2) then
              QQ_ALN_TEMP(ii,jj) = BLK%ULCT(i,j,k)
            else if (l==3) then
              QQ_ALN_TEMP(ii,jj) = BLK%VLCT(i,j,k)
            else if (l==4) then
              QQ_ALN_TEMP(ii,jj) = BLK%WLCT(i,j,k)
            else if (l==5) then
              QQ_ALN_TEMP(ii,jj) = BLK%PRSS(i,j,k)
            else
              QQ_ALN_TEMP(ii,jj) = BLK%YSPC(i,j,k,l-5)
            endif
          enddo
        enddo

        if (isc_reconstruction==1) then
          call weno5_align(QQ_ALN_TEMP,QL_ALN_TEMP,QR_ALN_TEMP,iimax,jjmax)
        else
          call muscl_align(QQ_ALN_TEMP,QL_ALN_TEMP,QR_ALN_TEMP,iorder_muscl,iimax,jjmax)
        endif

        do jj=1,jjmax+1
          do ii=1,iimax
            QL_ALN(ii,jj,l) = QL_ALN_TEMP(ii,jj)
            QR_ALN(ii,jj,l) = QR_ALN_TEMP(ii,jj)
          enddo
        enddo

      enddo

      deallocate(QQ_ALN_TEMP)
      deallocate(QL_ALN_TEMP)
      deallocate(QR_ALN_TEMP)

      call calc_roe

      do l=1,neqns
        do jj=1,kmax
          do ii=1,imax*jmax
            i = BLK%ista-1 + mod(ii-1,imax)+1
            j = BLK%jsta-1 + (ii-1)/imax+1
            k = BLK%ksta-1 + jj
            RHS(i,j,k,l) = RHS(i,j,k,l) - (FLUX_ALN(ii,nhalo1+jj+1,l)-FLUX_ALN(ii,nhalo1+jj,l))
          enddo
        enddo
      enddo

      deallocate(QL_ALN)
      deallocate(QR_ALN)
      deallocate(DXSX_ALN)
      deallocate(DXSY_ALN)
      deallocate(DXSZ_ALN)
      deallocate(XCRD_ALN)
      deallocate(YCRD_ALN)
      deallocate(ZCRD_ALN)
      deallocate(FLUX_ALN)
    end subroutine calc_zta
!-------------------------------------------------------------------------------
    subroutine calc_roe
      integer :: l,ii,jj,jl,jr
      real(8) :: sxm,sym,szm,sxyzm
      real(8) :: rhl,ul,vl,wl,pl,el,hl,akl,cl,ucl,hrt1l,hrt2l,hrt3l
      real(8) :: qql,ql1,ql2,ql3,ql4,ql5
      real(8) :: rhr,ur,vr,wr,pr,er,hr,akr,cr,ucr,hrt1r,hrt2r,hrt3r
      real(8) :: qqr,qr1,qr2,qr3,qr4,qr5
      real(8) :: rhm,um,vm,wm,hm,c2m,cm,akm,rwgt,wgtl,wgtr,ucm
      real(8) :: qqm,hrt1m,hrt2m,hrt3m
      real(8) :: ram1,ram2,ram3,eps1,eps2,eps3
      real(8) :: dq1,dq2,dq3,dq4,dq5
      real(8) :: qa1,qa2,qa3,qa4,qa5
      real(8) :: qb1,qb2,qb3,qb4,qb5
      real(8) :: alaa,albb,alcc,pbs,pbc,rsdu
      real(8) :: ysl,ysr
      real(8) :: eig1max,fl1,fr1,llfl1,llfr1,sig1,sigl1,sigr1,aal1,aar1

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
        jl = max0(1    ,jj-1)
        jr = min0(jjmax,jj  )
        do ii=1,iimax

          akl = cp_i/(cp_i-cr_i)

          rhl= QL_ALN(ii,jj,1)
          ul = QL_ALN(ii,jj,2)
          vl = QL_ALN(ii,jj,3)
          wl = QL_ALN(ii,jj,4)
          pl = QL_ALN(ii,jj,5)

          qql = 0.5d0*(ul*ul+vl*vl+wl*wl)

          el  = pl/(akl-1.0d0)+rhl*qql
          hl  = (el+pl)/rhl
          cl  = dsqrt((akl-1.0d0)*(hl-qql))

          ql1 = rhl
          ql2 = rhl*ul
          ql3 = rhl*vl
          ql4 = rhl*wl
          ql5 = el

          akr = cp_i/(cp_i-cr_i)

          rhr= QR_ALN(ii,jj,1)
          ur = QR_ALN(ii,jj,2)
          vr = QR_ALN(ii,jj,3)
          wr = QR_ALN(ii,jj,4)
          pr = QR_ALN(ii,jj,5)

          qqr = 0.5d0*(ur*ur+vr*vr+wr*wr)

          er = pr/(akr-1.0d0)+rhr*qqr
          hr = (er+pr)/rhr
          cr = dsqrt((akr-1.0d0)*(hr-qqr))

          qr1 = rhr
          qr2 = rhr*ur
          qr3 = rhr*vr
          qr4 = rhr*wr
          qr5 = er

          !Roe average
          sxm = DXSX_ALN(ii,jj)
          sym = DXSY_ALN(ii,jj)
          szm = DXSZ_ALN(ii,jj)
          sxyzm = dsqrt(sxm*sxm+sym*sym+szm*szm)

          rwgt= 1.0d0/(dsqrt(rhl)+dsqrt(rhr))
          wgtl= dsqrt(rhl)*rwgt
          wgtr= dsqrt(rhr)*rwgt
          rhm = dsqrt(rhl*rhr)
          um  = wgtl*ul+wgtr*ur
          vm  = wgtl*vl+wgtr*vr
          wm  = wgtl*wl+wgtr*wr
          hm  = wgtl*hl+wgtr*hr

          akm = 0.5d0*(akl+akr)
          c2m = (akm-1.0d0)*(hm-0.5d0*(um*um+vm*vm+wm*wm))
          cm  = dsqrt(c2m)

          qqm = 0.5d0*(um*um+vm*vm+wm*wm)

          !Rotation
          rxm = XCRD_ALN(ii,jj)-frame_origin(1)
          rym = YCRD_ALN(ii,jj)-frame_origin(2)
          rzm = ZCRD_ALN(ii,jj)-frame_origin(3)

          am = oy*rzm-oz*rym
          bm = oz*rxm-ox*rzm
          gm = ox*rym-oy*rxm
          acm= am*sxm+bm*sym+gm*szm

          ucl   = sxm*ul+sym*vl+szm*wl
          hrt1l = ucl-acm
          hrt2l = ucl-acm+cl*sxyzm
          hrt3l = ucl-acm-cl*sxyzm

          ucr   = sxm*ur+sym*vr+szm*wr
          hrt1r = ucr-acm
          hrt2r = ucr-acm+cr*sxyzm
          hrt3r = ucr-acm-cr*sxyzm

          ucm   = sxm*um+sym*vm+szm*wm
          hrt1m = ucm-acm
          hrt2m = ucm-acm+cm*sxyzm
          hrt3m = ucm-acm-cm*sxyzm

          ram1 = dabs(hrt1m)
          ram2 = dabs(hrt2m)
          ram3 = dabs(hrt3m)

          eps1 = dmax1(0.0d0,hrt1m-hrt1l,hrt1r-hrt1m)
          eps2 = dmax1(0.0d0,hrt2m-hrt2l,hrt2r-hrt2m)
          eps3 = dmax1(0.0d0,hrt3m-hrt3l,hrt3r-hrt3m)

          if (ram1.lt.eps1) ram1 = 0.5d0*(ram1*ram1/eps1+eps1)
          if (ram2.lt.eps2) ram2 = 0.5d0*(ram2*ram2/eps2+eps2)
          if (ram3.lt.eps3) ram3 = 0.5d0*(ram3*ram3/eps3+eps3)

          alaa = 0.5d0*(ram2-ram3)
          albb = 0.5d0*(ram2+ram3)-ram1
          alcc = ram1

          dq1 = qr1-ql1
          dq2 = qr2-ql2
          dq3 = qr3-ql3
          dq4 = qr4-ql4
          dq5 = qr5-ql5

          qa1 = 0.0d0
          qa2 = sxm
          qa3 = sym
          qa4 = szm
          qa5 = ucm

          qb1 = 1.0d0
          qb2 = um
          qb3 = vm
          qb4 = wm
          qb5 = hm

          pbs = (akm-1.0d0)*(dq1*qqm-(um*dq2+vm*dq3+wm*dq4)+dq5)
          pbc = pbs/(cm*cm)
          rsdu= sxm*dq2+sym*dq3+szm*dq4-dq1*ucm

          FLUX_ALN(ii,jj,1) = 0.5d0*(ql1*(ucl-acm)          +qr1*(ucr-acm)           &
          -alcc*dq1-alaa/(cm*sxyzm)*(pbs*qa1+rsdu*qb1)-albb*(pbc*qb1+rsdu/(sxyzm*sxyzm)*qa1))
          FLUX_ALN(ii,jj,2) = 0.5d0*(ql2*(ucl-acm)+pl*sxm   +qr2*(ucr-acm)+pr*sxm    &
          -alcc*dq2-alaa/(cm*sxyzm)*(pbs*qa2+rsdu*qb2)-albb*(pbc*qb2+rsdu/(sxyzm*sxyzm)*qa2))
          FLUX_ALN(ii,jj,3) = 0.5d0*(ql3*(ucl-acm)+pl*sym   +qr3*(ucr-acm)+pr*sym    &
          -alcc*dq3-alaa/(cm*sxyzm)*(pbs*qa3+rsdu*qb3)-albb*(pbc*qb3+rsdu/(sxyzm*sxyzm)*qa3))
          FLUX_ALN(ii,jj,4) = 0.5d0*(ql4*(ucl-acm)+pl*szm   +qr4*(ucr-acm)+pr*szm    &
          -alcc*dq4-alaa/(cm*sxyzm)*(pbs*qa4+rsdu*qb4)-albb*(pbc*qb4+rsdu/(sxyzm*sxyzm)*qa4))
          FLUX_ALN(ii,jj,5) = 0.5d0*(rhl*(ucl-acm)*hl+pl*acm+rhr*(ucr-acm)*hr+pr*acm &
          -alcc*dq5-alaa/(cm*sxyzm)*(pbs*qa5+rsdu*qb5)-albb*(pbc*qb5+rsdu/(sxyzm*sxyzm)*qa5))
        enddo
      enddo


      do l=1,neqns-5
        do jj=1,jjmax+1
          do ii=1,iimax

            rhl= QL_ALN(ii,jj,1)
            ul = QL_ALN(ii,jj,2)
            vl = QL_ALN(ii,jj,3)
            wl = QL_ALN(ii,jj,4)
            ysl= QL_ALN(ii,jj,5+l)

            rhr= QR_ALN(ii,jj,1)
            ur = QR_ALN(ii,jj,2)
            vr = QR_ALN(ii,jj,3)
            wr = QR_ALN(ii,jj,4)
            ysr= QR_ALN(ii,jj,5+l)

            sxm = DXSX_ALN(ii,jj)
            sym = DXSY_ALN(ii,jj)
            szm = DXSZ_ALN(ii,jj)

            !Rotation
            rxm = XCRD_ALN(ii,jj)-frame_origin(1)
            rym = YCRD_ALN(ii,jj)-frame_origin(2)
            rzm = ZCRD_ALN(ii,jj)-frame_origin(3)

            am = oy*rzm-oz*rym
            bm = oz*rxm-ox*rzm
            gm = ox*rym-oy*rxm
            acm= am*sxm+bm*sym+gm*szm

            ucl = sxm*ul+sym*vl+szm*wl
            ucr = sxm*ur+sym*vr+szm*wr

            eig1max = dmax1(dabs(ucl-acm),dabs(ucr-acm))

            ql1 = rhl*ysl
            fl1 = rhl*ysl*(ucl-acm)

            qr1 = rhr*ysr
            fr1 = rhr*ysr*(ucr-acm)

            llfl1 = 0.5d0*(fl1+eig1max*ql1)
            llfr1 = 0.5d0*(fr1-eig1max*qr1)

            sig1  = dsign(1.0d0,ucl-acm)*dsign(1.0d0,ucr-acm)
            sigl1 = 0.5d0*(1.0d0+dsign(1.0d0,ucl-acm))
            sigr1 = 0.5d0*(1.0d0-dsign(1.0d0,ucr-acm))

            aal1 = 0.5d0*((1.0d0+sig1)*sigl1*fl1+(1.0d0-sig1)*llfl1)
            aar1 = 0.5d0*((1.0d0+sig1)*sigr1*fr1+(1.0d0-sig1)*llfr1)

            FLUX_ALN(ii,jj,5+l) = aal1+aar1
          enddo
        enddo
      enddo

    end subroutine calc_roe

  end subroutine flux_riemann_roe

end module mdl_flux_riemann_roe
