module mdl_band
  implicit none

contains

  subroutine tridiag1(aa,bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8) :: bi,rbb
    real(8), dimension(1,jmax) :: tcc

    !Forward sweep
    rbb = 1.0d0/bb(1,1)
    tcc(1,1) = cc(1,1)*rbb
    do i=1,imax
      yy(i,1) = xx(i,1)*rbb
    enddo
    do j=1,jmax-1
      bi = 1.0d0/(bb(1,j+1)-aa(1,j+1)*tcc(1,j))
      tcc(1,j+1) = cc(1,j+1)*bi
      do i=1,imax
        yy(i,j+1) = (xx(i,j+1)-aa(1,j+1)*yy(i,j))*bi
      enddo
    enddo

    !Backward sweep
    do j=jmax,2,-1
      do i=1,imax
        yy(i,j-1) = yy(i,j-1)-tcc(1,j-1)*yy(i,j)
      enddo
    enddo

  end subroutine tridiag1

  subroutine tridiag2(aa,bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8) :: bi,rbb
    real(8), dimension(imax,jmax) :: tcc

    !Forward sweep
    do i=1,imax
      rbb = 1.0d0/bb(i,1)
      tcc(i,1) = cc(i,1)*rbb
      yy(i,1) = xx(i,1)*rbb
    enddo
    do j=1,jmax-1
      do i=1,imax
        bi = 1.0d0/(bb(i,j+1)-aa(i,j+1)*tcc(i,j))
        tcc(i,j+1) = cc(i,j+1)*bi
        yy(i,j+1) = (xx(i,j+1)-aa(i,j+1)*yy(i,j))*bi
      enddo
    enddo

    !Backward sweep
    do j=jmax,2,-1
      do i=1,imax
        yy(i,j-1) = yy(i,j-1)-tcc(i,j-1)*yy(i,j)
      enddo
    enddo

  end subroutine tridiag2

  subroutine tridiag1_cyc(aa,bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(1   ,jmax) :: at,bt,ct
    real(8), dimension(1   ,jmax) :: uu,qq
    real(8), dimension(imax,jmax) :: ww
    real(8) :: vn,rvq,vqq

    do j=1,jmax
      at(1,j) = aa(1,j)
      bt(1,j) = bb(1,j)
      ct(1,j) = cc(1,j)
      uu(1,j) = 0.0d0
    enddo
    bt(1,1) = 2.0d0*bb(1,1)
    bt(1,jmax) = bb(1,jmax)+cc(1,jmax)*aa(1,1)/bb(1,1)
    uu(1,1) = -bb(1,1)
    uu(1,jmax) = cc(1,jmax)

    call tridiag1(at,bt,ct,uu,qq,1,jmax)
    call tridiag1(at,bt,ct,xx,ww,imax,jmax)

    vn  = -aa(1,1)/bb(1,1)
    rvq = 1.0d0/(1.0d0+qq(1,1)+vn*qq(1,jmax))
    do j=1,jmax
      vqq = rvq*qq(1,j)
      do i=1,imax
        yy(i,j) = ww(i,j)-(ww(i,1)+vn*ww(i,jmax))*vqq
      enddo
    enddo

  end subroutine tridiag1_cyc

  subroutine tridiag2_cyc(aa,bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(imax,jmax) :: at,bt,ct
    real(8), dimension(imax,jmax) :: uu,qq
    real(8), dimension(imax,jmax) :: ww
    real(8) :: vn

    do j=1,jmax
      do i=1,imax
        at(i,j) = aa(i,j)
        bt(i,j) = bb(i,j)
        ct(i,j) = cc(i,j)
        uu(i,j) = 0.0d0
      enddo
    enddo

    do i=1,imax
      bt(i,1) = 2.0d0*bb(i,1)
      bt(i,jmax) = bb(i,jmax)+cc(i,jmax)*aa(i,1)/bb(i,1)
      uu(i,1) = -bb(i,1)
      uu(i,jmax) = cc(i,jmax)
    enddo

    call tridiag2(at,bt,ct,uu,qq,imax,jmax)
    call tridiag2(at,bt,ct,xx,ww,imax,jmax)

    do j=1,jmax
      do i=1,imax
        vn = -aa(i,1)/bb(i,1)
        yy(i,j) = ww(i,j) &
                 -(ww(i,1)+vn*ww(i,jmax))/(1.0d0+qq(i,1)+vn*qq(i,jmax))*qq(i,j)
      enddo
    enddo

  end subroutine tridiag2_cyc

  subroutine tridiag1_lower(aa,bb,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8) :: rbb

    !Forward sweep
    rbb = 1.0d0/bb(1,1)
    do i=1,imax
      yy(i,1) = xx(i,1)*rbb
    enddo

    do j=1,jmax-1
      rbb = 1.0d0/bb(1,j+1)
      do i=1,imax
        yy(i,j+1) = (xx(i,j+1)-aa(1,j+1)*yy(i,j))*rbb
      enddo
    enddo

  end subroutine tridiag1_lower

  subroutine tridiag2_lower(aa,bb,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j

    !Forward sweep
    do i=1,imax
      yy(i,1) = xx(i,1)/bb(i,1)
    enddo

    do j=1,jmax-1
      do i=1,imax
        yy(i,j+1) = (xx(i,j+1)-aa(i,j+1)*yy(i,j))/bb(i,j+1)
      enddo
    enddo

  end subroutine tridiag2_lower

  subroutine tridiag1_lower_cyc(aa,bb,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(1   ,jmax) :: at,bt
    real(8), dimension(1   ,jmax) :: uu,qq
    real(8), dimension(imax,jmax) :: ww
    real(8) :: vn,rvq,vqq

    do j=1,jmax
      at(1,j) = aa(1,j)
      bt(1,j) = bb(1,j)
      uu(1,j) = 0.0d0
    enddo
    bt(1,1) = 2.0d0*bb(1,1)
    uu(1,1) = -bb(1,1)

    call tridiag1_lower(at,bt,uu,qq,1   ,jmax)
    call tridiag1_lower(at,bt,xx,ww,imax,jmax)

    vn  = -aa(1,1)/bb(1,1)
    rvq = 1.0d0/(1.0d0+qq(1,1)+vn*qq(1,jmax))
    do j=1,jmax
      vqq = rvq*qq(1,j)
      do i=1,imax
        yy(i,j) = ww(i,j)-(ww(i,1)+vn*ww(i,jmax))*vqq
      enddo
    enddo

  end subroutine tridiag1_lower_cyc

  subroutine tridiag2_lower_cyc(aa,bb,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(imax,jmax) :: at,bt
    real(8), dimension(imax,jmax) :: uu,qq
    real(8), dimension(imax,jmax) :: ww
    real(8) :: vn

    do j=1,jmax
      do i=1,imax
        at(i,j) = aa(i,j)
        bt(i,j) = bb(i,j)
        uu(i,j) = 0.0d0
      enddo
    enddo
    do i=1,imax
      bt(i,1) = 2.0d0*bb(i,1)
      uu(i,1) = -bb(i,1)
    enddo

    call tridiag2_lower(at,bt,uu,qq,imax,jmax)
    call tridiag2_lower(at,bt,xx,ww,imax,jmax)

    do j=1,jmax
      do i=1,imax
        vn = -aa(i,1)/bb(i,1)
        yy(i,j) = ww(i,j) &
                 -(ww(i,1)+vn*ww(i,jmax))/(1.0d0+qq(i,1)+vn*qq(i,jmax))*qq(i,j)
      enddo
    enddo

  end subroutine tridiag2_lower_cyc

  subroutine tridiag1_upper(bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: bb,cc
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8) :: rbb

    !Backward sweep
    rbb = 1.0d0/bb(1,jmax)
    do i=1,imax
      yy(i,jmax) = xx(i,jmax)*rbb
    enddo

    do j=jmax,2,-1
      rbb = 1.0d0/bb(1,j-1)
      do i=1,imax
        yy(i,j-1) = (xx(i,j-1)-cc(1,j-1)*yy(i,j))*rbb
      enddo
    enddo

  end subroutine tridiag1_upper

  subroutine tridiag2_upper(bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: bb,cc,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j

    !Backward sweep
    do i=1,imax
      yy(i,jmax) = xx(i,jmax)/bb(i,jmax)
    enddo

    do j=jmax,2,-1
      do i=1,imax
        yy(i,j-1) = (xx(i,j-1)-cc(i,j-1)*yy(i,j))/bb(i,j-1)
      enddo
    enddo

  end subroutine tridiag2_upper

  subroutine tridiag1_upper_cyc(bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: bb,cc
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(1   ,jmax) :: bt,ct
    real(8), dimension(1   ,jmax) :: uu,qq
    real(8), dimension(imax,jmax) :: ww
    real(8) :: rvq,vqq

    do j=1,jmax
      bt(1,j) = bb(1,j)
      ct(1,j) = cc(1,j)
      uu(1,j) = 0.0d0
    enddo
    bt(1,1) = 2.0d0*bb(1,1)
    uu(1,1) = -bb(1,1)
    uu(1,jmax) = cc(1,jmax)

    call tridiag1_upper(bt,ct,uu,qq,1   ,jmax)
    call tridiag1_upper(bt,ct,xx,ww,imax,jmax)

    rvq = 1.0d0/(1.0d0+qq(1,1))
    do j=1,jmax
      vqq = rvq*qq(1,j)
      do i=1,imax
        yy(i,j) = ww(i,j)-ww(i,1)*vqq
      enddo
    enddo

  end subroutine tridiag1_upper_cyc

  subroutine tridiag2_upper_cyc(bb,cc,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: bb,cc,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(imax,jmax) :: bt,ct
    real(8), dimension(imax,jmax) :: uu,qq
    real(8), dimension(imax,jmax) :: ww

    do j=1,jmax
      do i=1,imax
        bt(i,j) = bb(i,j)
        ct(i,j) = cc(i,j)
        uu(i,j) = 0.0d0
      enddo
    enddo
    do i=1,imax
      bt(i,1) = 2.0d0*bb(i,1)
      uu(i,1) = -bb(i,1)
      uu(i,jmax) = cc(i,jmax)
    enddo

    call tridiag2_upper(bt,ct,uu,qq,imax,jmax)
    call tridiag2_upper(bt,ct,xx,ww,imax,jmax)

    do j=1,jmax
      do i=1,imax
        yy(i,j) = ww(i,j)-ww(i,1)/(1.0d0+qq(i,1))*qq(i,j)
      enddo
    enddo

  end subroutine tridiag2_upper_cyc

  subroutine pentadiag1(aa,bb,cc,dd,ee,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc,dd,ee
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8) :: ci,tbb,tcc
    real(8), dimension(1,jmax) :: tdd,tee

    tdd(1,1) = dd(1,1)/cc(1,1)
    tee(1,1) = ee(1,1)/cc(1,1)
    do i=1,imax
      yy(i,1) = xx(i,1)/cc(1,1)
      yy(i,2) = xx(i,2)
    enddo
    tbb = bb(1,2)
    tcc = cc(1,2)

    do j=1,jmax-2
      ci = 1.0d0/(tcc-tbb*tdd(1,j))
      tdd(1,j+1) = (dd(1,j+1)-tbb*tee(1,j))*ci
      tee(1,j+1) = ee(1,j+1)*ci
      do i=1,imax
        yy(i,j+1) = (yy(i,j+1)-tbb*yy(i,j))*ci
      enddo

      tbb = bb(1,j+2)-aa(1,j+2)*tdd(1,j)
      tcc = cc(1,j+2)-aa(1,j+2)*tee(1,j)
      do i=1,imax
        yy(i,j+2) = xx(i,j+2)-aa(1,j+2)*yy(i,j)
      enddo
    enddo

    ci = 1.0d0/(tcc-tbb*tdd(1,jmax-1))
    do i=1,imax
      yy(i,jmax) = (yy(i,jmax)-tbb*yy(i,jmax-1))*ci
    enddo

    do j=jmax,3,-1
      do i=1,imax
        yy(i,j-1) = yy(i,j-1)-tdd(1,j-1)*yy(i,j)
        yy(i,j-2) = yy(i,j-2)-tee(1,j-2)*yy(i,j)
      enddo
    enddo

    do i=1,imax
      yy(i,1) = yy(i,1)-tdd(1,1)*yy(i,2)
    enddo

  end subroutine pentadiag1

  subroutine pentadiag2(aa,bb,cc,dd,ee,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc,dd,ee,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8) :: ci
    real(8), dimension(imax)      :: tbb,tcc
    real(8), dimension(imax,jmax) :: tdd,tee

    do i=1,imax
      tdd(i,1) = dd(i,1)/cc(i,1)
      tee(i,1) = ee(i,1)/cc(i,1)
      yy(i,1) = xx(i,1)/cc(i,1)
      yy(i,2) = xx(i,2)

      tbb(i) = bb(i,2)
      tcc(i) = cc(i,2)
    enddo

    do j=1,jmax-2
      do i=1,imax
        ci = 1.0d0/(tcc(i)-tbb(i)*tdd(i,j))
        tdd(i,j+1) = (dd(i,j+1)-tbb(i)*tee(i,j))*ci
        tee(i,j+1) = ee(i,j+1)*ci
        yy(i,j+1) = (yy(i,j+1)-tbb(i)*yy(i,j))*ci

        tbb(i) = bb(i,j+2)-aa(i,j+2)*tdd(i,j)
        tcc(i) = cc(i,j+2)-aa(i,j+2)*tee(i,j)
        yy(i,j+2) = xx(i,j+2)-aa(i,j+2)*yy(i,j)
      enddo
    enddo

    do i=1,imax
      ci = 1.0d0/(tcc(i)-tbb(i)*tdd(i,jmax-1))
      yy(i,jmax) = (yy(i,jmax)-tbb(i)*yy(i,jmax-1))*ci
    enddo

    do j=jmax,3,-1
      do i=1,imax
        yy(i,j-1) = yy(i,j-1)-tdd(i,j-1)*yy(i,j)
        yy(i,j-2) = yy(i,j-2)-tee(i,j-2)*yy(i,j)
      enddo
    enddo

    do i=1,imax
      yy(i,1) = yy(i,1)-tdd(i,1)*yy(i,2)
    enddo

  end subroutine pentadiag2

  subroutine pentadiag1_cyc(aa,bb,cc,dd,ee,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc,dd,ee
    real(8), dimension(:,:), intent(in ) :: xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(1   ,jmax) :: u1,u2,u3,u4
    real(8), dimension(1   ,jmax) :: z1,z2,z3,z4
    real(8), dimension(imax,jmax) :: ww
    real(8), dimension(1,4,4)     :: p,q
    real(8) :: r1,r2,r3,r4
    real(8), dimension(imax) :: s1,s2,s3,s4

    do j=1,jmax
      u1(1,j) = 0.0d0
      u2(1,j) = 0.0d0
      u3(1,j) = 0.0d0
      u4(1,j) = 0.0d0
    enddo

    u1(1,1) = 1.0d0
    u2(1,2) = 1.0d0
    u3(1,jmax-1) = 1.0d0
    u4(1,jmax) = 1.0d0

    call pentadiag1(aa,bb,cc,dd,ee,u1,z1,1,jmax)
    call pentadiag1(aa,bb,cc,dd,ee,u2,z2,1,jmax)
    call pentadiag1(aa,bb,cc,dd,ee,u3,z3,1,jmax)
    call pentadiag1(aa,bb,cc,dd,ee,u4,z4,1,jmax)
    call pentadiag1(aa,bb,cc,dd,ee,xx,ww,imax,jmax)

    p(1,1,1) = aa(1,1)*z1(1,jmax-1) + bb(1,1)*z1(1,jmax) + 1.0d0
    p(1,1,2) = aa(1,1)*z2(1,jmax-1) + bb(1,1)*z2(1,jmax)
    p(1,1,3) = aa(1,1)*z3(1,jmax-1) + bb(1,1)*z3(1,jmax)
    p(1,1,4) = aa(1,1)*z4(1,jmax-1) + bb(1,1)*z4(1,jmax)
    !---------------------------------------------------
    p(1,2,1) = aa(1,2)*z1(1,jmax)
    p(1,2,2) = aa(1,2)*z2(1,jmax) + 1.0d0
    p(1,2,3) = aa(1,2)*z3(1,jmax)
    p(1,2,4) = aa(1,2)*z4(1,jmax)
    !---------------------------------------------------
    p(1,3,1) = ee(1,jmax-1)*z1(1,1)
    p(1,3,2) = ee(1,jmax-1)*z2(1,1)
    p(1,3,3) = ee(1,jmax-1)*z3(1,1) + 1.0d0
    p(1,3,4) = ee(1,jmax-1)*z4(1,1)
    !---------------------------------------------------
    p(1,4,1) = dd(1,jmax)*z1(1,1) + ee(1,jmax)*z1(1,2)
    p(1,4,2) = dd(1,jmax)*z2(1,1) + ee(1,jmax)*z2(1,2)
    p(1,4,3) = dd(1,jmax)*z3(1,1) + ee(1,jmax)*z3(1,2)
    p(1,4,4) = dd(1,jmax)*z4(1,1) + ee(1,jmax)*z4(1,2) + 1.0d0
    !---------------------------------------------------

    call inv44(p,q,1)

    do i=1,imax
      r1 = aa(1,1)*ww(i,jmax-1) + bb(1,1)*ww(i,jmax)
      r2 = aa(1,2)*ww(i,jmax)
      r3 = ee(1,jmax-1)*ww(i,1)
      r4 = dd(1,jmax)*ww(i,1) + ee(1,jmax)*ww(i,2)

      s1(i) = q(1,1,1)*r1+q(1,1,2)*r2+q(1,1,3)*r3+q(1,1,4)*r4
      s2(i) = q(1,2,1)*r1+q(1,2,2)*r2+q(1,2,3)*r3+q(1,2,4)*r4
      s3(i) = q(1,3,1)*r1+q(1,3,2)*r2+q(1,3,3)*r3+q(1,3,4)*r4
      s4(i) = q(1,4,1)*r1+q(1,4,2)*r2+q(1,4,3)*r3+q(1,4,4)*r4
    enddo

    do j=1,jmax
      do i=1,imax
        yy(i,j) = ww(i,j) &
                 -(z1(1,j)*s1(i)+z2(1,j)*s2(i)+z3(1,j)*s3(i)+z4(1,j)*s4(i))
      enddo
    enddo

  end subroutine pentadiag1_cyc

  subroutine pentadiag2_cyc(aa,bb,cc,dd,ee,xx,yy,imax,jmax)
    integer, intent(in) :: imax,jmax
    real(8), dimension(:,:), intent(in ) :: aa,bb,cc,dd,ee,xx
    real(8), dimension(:,:), intent(out) :: yy

    integer :: i,j
    real(8), dimension(imax,jmax) :: u1,u2,u3,u4
    real(8), dimension(imax,jmax) :: z1,z2,z3,z4
    real(8), dimension(imax,jmax) :: ww
    real(8), dimension(imax,4,4)  :: p,q
    real(8) :: r1,r2,r3,r4
    real(8), dimension(imax) :: s1,s2,s3,s4

    do j=1,jmax
      do i=1,imax
        u1(i,j) = 0.0d0
        u2(i,j) = 0.0d0
        u3(i,j) = 0.0d0
        u4(i,j) = 0.0d0
      enddo
    enddo

    do i=1,imax
      u1(i,1) = 1.0d0
      u2(i,2) = 1.0d0
      u3(i,jmax-1) = 1.0d0
      u4(i,jmax) = 1.0d0
    enddo

    call pentadiag2(aa,bb,cc,dd,ee,u1,z1,imax,jmax)
    call pentadiag2(aa,bb,cc,dd,ee,u2,z2,imax,jmax)
    call pentadiag2(aa,bb,cc,dd,ee,u3,z3,imax,jmax)
    call pentadiag2(aa,bb,cc,dd,ee,u4,z4,imax,jmax)
    call pentadiag2(aa,bb,cc,dd,ee,xx,ww,imax,jmax)

    do i=1,imax
      p(i,1,1) = aa(i,1)*z1(i,jmax-1) + bb(i,1)*z1(i,jmax) + 1.0d0
      p(i,1,2) = aa(i,1)*z2(i,jmax-1) + bb(i,1)*z2(i,jmax)
      p(i,1,3) = aa(i,1)*z3(i,jmax-1) + bb(i,1)*z3(i,jmax)
      p(i,1,4) = aa(i,1)*z4(i,jmax-1) + bb(i,1)*z4(i,jmax)
      !---------------------------------------------------
      p(i,2,1) = aa(i,2)*z1(i,jmax)
      p(i,2,2) = aa(i,2)*z2(i,jmax) + 1.0d0
      p(i,2,3) = aa(i,2)*z3(i,jmax)
      p(i,2,4) = aa(i,2)*z4(i,jmax)
      !---------------------------------------------------
      p(i,3,1) = ee(i,jmax-1)*z1(i,1)
      p(i,3,2) = ee(i,jmax-1)*z2(i,1)
      p(i,3,3) = ee(i,jmax-1)*z3(i,1) + 1.0d0
      p(i,3,4) = ee(i,jmax-1)*z4(i,1)
      !---------------------------------------------------
      p(i,4,1) = dd(i,jmax)*z1(i,1) + ee(i,jmax)*z1(i,2)
      p(i,4,2) = dd(i,jmax)*z2(i,1) + ee(i,jmax)*z2(i,2)
      p(i,4,3) = dd(i,jmax)*z3(i,1) + ee(i,jmax)*z3(i,2)
      p(i,4,4) = dd(i,jmax)*z4(i,1) + ee(i,jmax)*z4(i,2) + 1.0d0
    enddo

    call inv44(p,q,imax)

    do i=1,imax
      r1 = aa(i,1)*ww(i,jmax-1) + bb(i,1)*ww(i,jmax)
      r2 = aa(i,2)*ww(i,jmax)
      r3 = ee(i,jmax-1)*ww(i,1)
      r4 = dd(i,jmax)*ww(i,1) + ee(i,jmax)*ww(i,2)

      s1(i) = q(i,1,1)*r1+q(i,1,2)*r2+q(i,1,3)*r3+q(i,1,4)*r4
      s2(i) = q(i,2,1)*r1+q(i,2,2)*r2+q(i,2,3)*r3+q(i,2,4)*r4
      s3(i) = q(i,3,1)*r1+q(i,3,2)*r2+q(i,3,3)*r3+q(i,3,4)*r4
      s4(i) = q(i,4,1)*r1+q(i,4,2)*r2+q(i,4,3)*r3+q(i,4,4)*r4
    enddo

    do j=1,jmax
      do i=1,imax
        yy(i,j) = ww(i,j) &
                 -(z1(i,j)*s1(i)+z2(i,j)*s2(i)+z3(i,j)*s3(i)+z4(i,j)*s4(i))
      enddo
    enddo

  end subroutine pentadiag2_cyc

  subroutine inv44(pp,qq,imax)
    integer, intent(in) :: imax
    real(8), dimension(:,:,:), intent(in ) :: pp
    real(8), dimension(:,:,:), intent(out) :: qq

    integer :: i
    real(8) :: rpp
    real(8), dimension(imax,4,4) :: rr

    do i=1,imax
      rpp = 1.0d0/(pp(i,1,1))
      rr(i,1,2) = pp(i,1,2)*rpp
      rr(i,1,3) = pp(i,1,3)*rpp
      rr(i,1,4) = pp(i,1,4)*rpp
      qq(i,1,1) = rpp

      rr(i,2,2) = pp(i,2,2)-pp(i,2,1)*rr(i,1,2)
      rr(i,2,3) = pp(i,2,3)-pp(i,2,1)*rr(i,1,3)
      rr(i,2,4) = pp(i,2,4)-pp(i,2,1)*rr(i,1,4)
      qq(i,2,1) =-pp(i,2,1)*qq(i,1,1)

      rr(i,3,2) = pp(i,3,2)-pp(i,3,1)*rr(i,1,2)
      rr(i,3,3) = pp(i,3,3)-pp(i,3,1)*rr(i,1,3)
      rr(i,3,4) = pp(i,3,4)-pp(i,3,1)*rr(i,1,4)
      qq(i,3,1) =-pp(i,3,1)*qq(i,1,1)

      rr(i,4,2) = pp(i,4,2)-pp(i,4,1)*rr(i,1,2)
      rr(i,4,3) = pp(i,4,3)-pp(i,4,1)*rr(i,1,3)
      rr(i,4,4) = pp(i,4,4)-pp(i,4,1)*rr(i,1,4)
      qq(i,4,1) =-pp(i,4,1)*qq(i,1,1)
      !-------------------------------------------------
      rpp = 1.0d0/(rr(i,2,2))
      rr(i,2,3) = rr(i,2,3)*rpp
      rr(i,2,4) = rr(i,2,4)*rpp
      qq(i,2,1) = qq(i,2,1)*rpp
      qq(i,2,2) = rpp

      rr(i,3,3) = rr(i,3,3)-rr(i,3,2)*rr(i,2,3)
      rr(i,3,4) = rr(i,3,4)-rr(i,3,2)*rr(i,2,4)
      qq(i,3,1) = qq(i,3,1)-rr(i,3,2)*qq(i,2,1)
      qq(i,3,2) =-rr(i,3,2)*qq(i,2,2)

      rr(i,4,3) = rr(i,4,3)-rr(i,4,2)*rr(i,2,3)
      rr(i,4,4) = rr(i,4,4)-rr(i,4,2)*rr(i,2,4)
      qq(i,4,1) = qq(i,4,1)-rr(i,4,2)*qq(i,2,1)
      qq(i,4,2) =-rr(i,4,2)*qq(i,2,2)
      !-------------------------------------------------
      rpp = 1.0d0/(rr(i,3,3))
      rr(i,3,4) = rr(i,3,4)*rpp
      qq(i,3,1) = qq(i,3,1)*rpp
      qq(i,3,2) = qq(i,3,2)*rpp
      qq(i,3,3) = rpp

      rr(i,4,4) = rr(i,4,4)-rr(i,4,3)*rr(i,3,4)
      qq(i,4,1) = qq(i,4,1)-rr(i,4,3)*qq(i,3,1)
      qq(i,4,2) = qq(i,4,2)-rr(i,4,3)*qq(i,3,2)
      qq(i,4,3) =-rr(i,4,3)*qq(i,3,3)
      !-------------------------------------------------
      rpp = 1.0d0/(rr(i,4,4))
      qq(i,4,1) = qq(i,4,1)*rpp
      qq(i,4,2) = qq(i,4,2)*rpp
      qq(i,4,3) = qq(i,4,3)*rpp
      qq(i,4,4) = rpp
      !-------------------------------------------------
      qq(i,3,1) = qq(i,3,1)-rr(i,3,4)*qq(i,4,1)
      qq(i,3,2) = qq(i,3,2)-rr(i,3,4)*qq(i,4,2)
      qq(i,3,3) = qq(i,3,3)-rr(i,3,4)*qq(i,4,3)
      qq(i,3,4) =-rr(i,3,4)*qq(i,4,4)

      qq(i,2,1) = qq(i,2,1)-rr(i,2,4)*qq(i,4,1)
      qq(i,2,2) = qq(i,2,2)-rr(i,2,4)*qq(i,4,2)
      qq(i,2,3) =-rr(i,2,4)*qq(i,4,3)
      qq(i,2,4) =-rr(i,2,4)*qq(i,4,4)

      qq(i,1,1) = qq(i,1,1)-rr(i,1,4)*qq(i,4,1)
      qq(i,1,2) =-rr(i,1,4)*qq(i,4,2)
      qq(i,1,3) =-rr(i,1,4)*qq(i,4,3)
      qq(i,1,4) =-rr(i,1,4)*qq(i,4,4)
      !-------------------------------------------------
      qq(i,2,1) = qq(i,2,1)-rr(i,2,3)*qq(i,3,1)
      qq(i,2,2) = qq(i,2,2)-rr(i,2,3)*qq(i,3,2)
      qq(i,2,3) = qq(i,2,3)-rr(i,2,3)*qq(i,3,3)
      qq(i,2,4) = qq(i,2,4)-rr(i,2,3)*qq(i,3,4)

      qq(i,1,1) = qq(i,1,1)-rr(i,1,3)*qq(i,3,1)
      qq(i,1,2) = qq(i,1,2)-rr(i,1,3)*qq(i,3,2)
      qq(i,1,3) = qq(i,1,3)-rr(i,1,3)*qq(i,3,3)
      qq(i,1,4) = qq(i,1,4)-rr(i,1,3)*qq(i,3,4)
      !-------------------------------------------------
      qq(i,1,1) = qq(i,1,1)-rr(i,1,2)*qq(i,2,1)
      qq(i,1,2) = qq(i,1,2)-rr(i,1,2)*qq(i,2,2)
      qq(i,1,3) = qq(i,1,3)-rr(i,1,2)*qq(i,2,3)
      qq(i,1,4) = qq(i,1,4)-rr(i,1,2)*qq(i,2,4)
    enddo

  end subroutine inv44

end module mdl_band
