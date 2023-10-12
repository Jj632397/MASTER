module mdl_lp_filter
  use mdl_band
  implicit none

  integer, private, parameter :: align_mem_access = 1 !<0>:simple <1>:longer vector length

contains
!-------------------------------------------------------------------------------
!idir
!   <1>:xsi <2>:eta <3>:zta
!
!iperi
!   <1>:periodic <else>:bounded
!   If 'iperi' is '1', isc(6) with afs(6) is applied for all points.
!
!icompact
!   <1>:implicit filter <else>:explicit filter
!   If expricit filter is used, afs(1)-afs(11) are zero.
!
! isc(1),isc(11)  [Point 1,N]
!   <any>:0th-order
! isc(2),isc(10)  [Point 2,N-1]
!   <2>:2nd-order <4>:4th-order <6>:6th-order <8>:8th-order <10>:10th-order
!   <else>:0th-order
! isc(3),isc(9)   [Point 3,N-2]
!   <4>:4th-order <6>:6th-order <8>:8th-order <10>:10th-order
!   <else>:0th-order
! isc(4),isc(8)   [Point 4,N-3]
!   <6>:6th-order <8>:8th-order <10>:10th-order
!   <else>:0th-order
! isc(5),isc(7)   [Point 5,N-4]
!   <8>:8th-order <10>:10th-order
!   <else>:0th-order
! isc(6)          [Interior]
!   <8>:8th-order <10>:10th-order
!   <else>:0th-order
!
!Free Parameter -0.5 < af < 0.5
! afs(1)  [Point 1]
! afs(2)  [Point 2]
! afs(3)  [Point 3]
! afs(4)  [Point 4]
! afs(5)  [Point 5]
! afs(6)  [Interior]
! afs(7)  [Point N-4]
! afs(8)  [Point N-3]
! afs(9)  [Point N-2]
! afs(10) [Point N-1]
! afs(11) [Point N]
!-------------------------------------------------------------------------------
  subroutine calc_lp_filter(F,idir,iperi,icompact,isc,afs,imax,jmax,kmax)
    real(8), dimension(:,:,:), intent(inout) :: F
    integer, intent(in) :: idir,iperi,icompact
    integer, dimension(:), intent(in) :: isc
    real(8), dimension(:), intent(in) :: afs
    integer, intent(in) :: imax,jmax,kmax

    integer :: i,j,k
    integer :: ii,jj,iimax,jjmax

    real(8), allocatable, dimension(:,:) :: F_ALN

    if (idir==1) then
      iimax= jmax*kmax
      jjmax= imax
    else if (idir==2) then
      iimax= imax*kmax
      jjmax= jmax
    else if (idir==3) then
      iimax= imax*jmax
      jjmax= kmax
    endif

    allocate(F_ALN(iimax,jjmax))

    if (idir==1) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = j+jmax*(k-1)
              jj = i
              F_ALN(ii,jj) = F(i,j,k)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = jj
            j = mod(ii-1,jmax)+1
            k = (ii-1)/jmax+1
            F_ALN(ii,jj) = F(i,j,k)
          enddo
        enddo
      endif

    else if (idir==2) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(k-1)
              jj = j
              F_ALN(ii,jj) = F(i,j,k)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = jj
            k = (ii-1)/imax+1
            F_ALN(ii,jj) = F(i,j,k)
          enddo
        enddo
      endif

    else if (idir==3) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(j-1)
              jj = k
              F_ALN(ii,jj) = F(i,j,k)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = (ii-1)/imax+1
            k = jj
            F_ALN(ii,jj) = F(i,j,k)
          enddo
        enddo
      endif

    endif


    call calc_lp_filter_align(F_ALN,iperi,icompact,isc,afs,iimax,jjmax)


    if (idir==1) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = j+jmax*(k-1)
              jj = i
              F(i,j,k) = F_ALN(ii,jj)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = jj
            j = mod(ii-1,jmax)+1
            k = (ii-1)/jmax+1
            F(i,j,k) = F_ALN(ii,jj)
          enddo
        enddo
      endif

    else if (idir==2) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(k-1)
              jj = j
              F(i,j,k) = F_ALN(ii,jj)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = jj
            k = (ii-1)/imax+1
            F(i,j,k) = F_ALN(ii,jj)
          enddo
        enddo
      endif

    else if (idir==3) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(j-1)
              jj = k
              F(i,j,k) = F_ALN(ii,jj)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = (ii-1)/imax+1
            k = jj
            F(i,j,k) = F_ALN(ii,jj)
          enddo
        enddo
      endif

    endif


    deallocate(F_ALN)

  end subroutine calc_lp_filter
!-------------------------------------------------------------------------------
  subroutine calc_lp_filter_align(F,iperi,icompact,isc,afs,iimax,jjmax)
    real(8), dimension(:,:), intent(inout) :: F
    integer, intent(in) :: iperi,icompact
    integer, dimension(:), intent(in) :: isc
    real(8), dimension(:), intent(in) :: afs
    integer, intent(in) :: iimax,jjmax

    integer :: ii,jj
    real(8), dimension(:,:), allocatable :: aa,bb,cc,xx
    real(8) :: af,a0,a1,a2,a3,a4,a5
    real(8) :: af_1 ,a1_1 ,a2_1 ,a3_1 ,a4_1 ,a5_1 ,a6_1 ,a7_1 ,a8_1 ,a9_1, a10_1 ,a11_1
    real(8) :: af_2 ,a1_2 ,a2_2 ,a3_2 ,a4_2 ,a5_2 ,a6_2 ,a7_2 ,a8_2 ,a9_2, a10_2 ,a11_2
    real(8) :: af_3 ,a1_3 ,a2_3 ,a3_3 ,a4_3 ,a5_3 ,a6_3 ,a7_3 ,a8_3 ,a9_3, a10_3 ,a11_3
    real(8) :: af_4 ,a1_4 ,a2_4 ,a3_4 ,a4_4 ,a5_4 ,a6_4 ,a7_4 ,a8_4 ,a9_4, a10_4 ,a11_4
    real(8) :: af_5 ,a1_5 ,a2_5 ,a3_5 ,a4_5 ,a5_5 ,a6_5 ,a7_5 ,a8_5 ,a9_5, a10_5 ,a11_5
    real(8) :: af_n ,a1_n ,a2_n ,a3_n ,a4_n ,a5_n ,a6_n ,a7_n ,a8_n ,a9_n, a10_n ,a11_n
    real(8) :: af_n1,a1_n1,a2_n1,a3_n1,a4_n1,a5_n1,a6_n1,a7_n1,a8_n1,a9_n1,a10_n1,a11_n1
    real(8) :: af_n2,a1_n2,a2_n2,a3_n2,a4_n2,a5_n2,a6_n2,a7_n2,a8_n2,a9_n2,a10_n2,a11_n2
    real(8) :: af_n3,a1_n3,a2_n3,a3_n3,a4_n3,a5_n3,a6_n3,a7_n3,a8_n3,a9_n3,a10_n3,a11_n3
    real(8) :: af_n4,a1_n4,a2_n4,a3_n4,a4_n4,a5_n4,a6_n4,a7_n4,a8_n4,a9_n4,a10_n4,a11_n4

    if (jjmax<11) return

    call set_coefficients

    allocate(aa(1    ,jjmax))
    allocate(bb(1    ,jjmax))
    allocate(cc(1    ,jjmax))
    allocate(xx(iimax,jjmax))

    !Set L.H.S
    if (iperi==1) then
      do jj=1,jjmax
        aa(1,jj) = af
        bb(1,jj) = 1.0d0
        cc(1,jj) = af
      enddo
    else
      !Point 1
      aa(1,1) = 0.0d0
      bb(1,1) = 1.0d0
      cc(1,1) = af_1
      !Point 2
      aa(1,2) = af_2
      bb(1,2) = 1.0d0
      cc(1,2) = af_2
      !Point 3
      aa(1,3) = af_3
      bb(1,3) = 1.0d0
      cc(1,3) = af_3
      !Point 4
      aa(1,4) = af_4
      bb(1,4) = 1.0d0
      cc(1,4) = af_4
      !Point 5
      aa(1,5) = af_5
      bb(1,5) = 1.0d0
      cc(1,5) = af_5
      !Point 6 - N-5
      do jj=6,jjmax-5
        aa(1,jj) = af
        bb(1,jj) = 1.0d0
        cc(1,jj) = af
      enddo
      !Point N-4
      aa(1,jjmax-4) = af_n4
      bb(1,jjmax-4) = 1.0d0
      cc(1,jjmax-4) = af_n4
      !Point N-3
      aa(1,jjmax-3) = af_n3
      bb(1,jjmax-3) = 1.0d0
      cc(1,jjmax-3) = af_n3
      !Point N-2
      aa(1,jjmax-2) = af_n2
      bb(1,jjmax-2) = 1.0d0
      cc(1,jjmax-2) = af_n2
      !Point N-1
      aa(1,jjmax-1) = af_n1
      bb(1,jjmax-1) = 1.0d0
      cc(1,jjmax-1) = af_n1
      !Point N
      aa(1,jjmax  ) = af_n
      bb(1,jjmax  ) = 1.0d0
      cc(1,jjmax  ) = 0.0d0
    endif

    call set_rhs

    if (icompact==1) then
      if (iperi==1) then
        call tridiag1_cyc(aa,bb,cc,xx,F,iimax,jjmax)
      else
        call tridiag1(aa,bb,cc,xx,F,iimax,jjmax)
      endif
    else
      F = xx
    endif

    deallocate(aa)
    deallocate(bb)
    deallocate(cc)
    deallocate(xx)

  contains
!-------------------------------------------------------------------------------
    subroutine set_rhs

      do jj=6,jjmax-5
        do ii=1,iimax
          xx(ii,jj) = a0* F(ii,jj  )             &
                     +a1*(F(ii,jj+1)+F(ii,jj-1)) &
                     +a2*(F(ii,jj+2)+F(ii,jj-2)) &
                     +a3*(F(ii,jj+3)+F(ii,jj-3)) &
                     +a4*(F(ii,jj+4)+F(ii,jj-4)) &
                     +a5*(F(ii,jj+5)+F(ii,jj-5))
        enddo
      enddo

      if (iperi==1) then

        do ii=1,iimax
          !Point 1
          xx(ii,1) = a0* F(ii,1)                &
                    +a1*(F(ii,2)+F(ii,jjmax  )) &
                    +a2*(F(ii,3)+F(ii,jjmax-1)) &
                    +a3*(F(ii,4)+F(ii,jjmax-2)) &
                    +a4*(F(ii,5)+F(ii,jjmax-3)) &
                    +a5*(F(ii,6)+F(ii,jjmax-4))
          !Point 2
          xx(ii,2) = a0* F(ii,2)                &
                    +a1*(F(ii,3)+F(ii,1      )) &
                    +a2*(F(ii,4)+F(ii,jjmax  )) &
                    +a3*(F(ii,5)+F(ii,jjmax-1)) &
                    +a4*(F(ii,6)+F(ii,jjmax-2)) &
                    +a5*(F(ii,7)+F(ii,jjmax-3))
          !Point 3
          xx(ii,3) = a0* F(ii,3)                &
                    +a1*(F(ii,4)+F(ii,2      )) &
                    +a2*(F(ii,5)+F(ii,1      )) &
                    +a3*(F(ii,6)+F(ii,jjmax  )) &
                    +a4*(F(ii,7)+F(ii,jjmax-1)) &
                    +a5*(F(ii,8)+F(ii,jjmax-2))
          !Point 4
          xx(ii,4) = a0* F(ii,4)                &
                    +a1*(F(ii,5)+F(ii,3      )) &
                    +a2*(F(ii,6)+F(ii,2      )) &
                    +a3*(F(ii,7)+F(ii,1      )) &
                    +a4*(F(ii,8)+F(ii,jjmax  )) &
                    +a5*(F(ii,9)+F(ii,jjmax-1))
          !Point 5
          xx(ii,5) = a0* F(ii,5 )              &
                    +a1*(F(ii,6 )+F(ii,4    )) &
                    +a2*(F(ii,7 )+F(ii,3    )) &
                    +a3*(F(ii,8 )+F(ii,2    )) &
                    +a4*(F(ii,9 )+F(ii,1    )) &
                    +a5*(F(ii,10)+F(ii,jjmax))


          !Point N
          xx(ii,jjmax  ) = a0* F(ii,jjmax)            &
                          +a1*(F(ii,1)+F(ii,jjmax-1)) &
                          +a2*(F(ii,2)+F(ii,jjmax-2)) &
                          +a3*(F(ii,3)+F(ii,jjmax-3)) &
                          +a4*(F(ii,4)+F(ii,jjmax-4)) &
                          +a5*(F(ii,5)+F(ii,jjmax-5))
          !Point N-1
          xx(ii,jjmax-1) = a0* F(ii,jjmax-1)              &
                          +a1*(F(ii,jjmax)+F(ii,jjmax-2)) &
                          +a2*(F(ii,1    )+F(ii,jjmax-3)) &
                          +a3*(F(ii,2    )+F(ii,jjmax-4)) &
                          +a4*(F(ii,3    )+F(ii,jjmax-5)) &
                          +a5*(F(ii,4    )+F(ii,jjmax-6))
          !Point N-2
          xx(ii,jjmax-2) = a0* F(ii,jjmax-2)                &
                          +a1*(F(ii,jjmax-1)+F(ii,jjmax-3)) &
                          +a2*(F(ii,jjmax  )+F(ii,jjmax-4)) &
                          +a3*(F(ii,1      )+F(ii,jjmax-5)) &
                          +a4*(F(ii,2      )+F(ii,jjmax-6)) &
                          +a5*(F(ii,3      )+F(ii,jjmax-7))
          !Point N-3
          xx(ii,jjmax-3) = a0* F(ii,jjmax-3)                &
                          +a1*(F(ii,jjmax-2)+F(ii,jjmax-4)) &
                          +a2*(F(ii,jjmax-1)+F(ii,jjmax-5)) &
                          +a3*(F(ii,jjmax  )+F(ii,jjmax-6)) &
                          +a4*(F(ii,1      )+F(ii,jjmax-7)) &
                          +a5*(F(ii,2      )+F(ii,jjmax-8))
          !Point N-4
          xx(ii,jjmax-4) = a0* F(ii,jjmax-4)                &
                          +a1*(F(ii,jjmax-3)+F(ii,jjmax-5)) &
                          +a2*(F(ii,jjmax-2)+F(ii,jjmax-6)) &
                          +a3*(F(ii,jjmax-1)+F(ii,jjmax-7)) &
                          +a4*(F(ii,jjmax  )+F(ii,jjmax-8)) &
                          +a5*(F(ii,1      )+F(ii,jjmax-9))
        enddo

      else

        do ii=1,iimax
          xx(ii,1) = a1_1 *F(ii,1 )+a2_1 *F(ii,2 ) &
                    +a3_1 *F(ii,3 )+a4_1 *F(ii,4 ) &
                    +a5_1 *F(ii,5 )+a6_1 *F(ii,6 ) &
                    +a7_1 *F(ii,7 )+a8_1 *F(ii,8 ) &
                    +a9_1 *F(ii,9 )+a10_1*F(ii,10) &
                    +a11_1*F(ii,11)
          xx(ii,2) = a1_2 *F(ii,1 )+a2_2 *F(ii,2 ) &
                    +a3_2 *F(ii,3 )+a4_2 *F(ii,4 ) &
                    +a5_2 *F(ii,5 )+a6_2 *F(ii,6 ) &
                    +a7_2 *F(ii,7 )+a8_2 *F(ii,8 ) &
                    +a9_2 *F(ii,9 )+a10_2*F(ii,10) &
                    +a11_2*F(ii,11)
          xx(ii,3) = a1_3 *F(ii,1 )+a2_3 *F(ii,2 ) &
                    +a3_3 *F(ii,3 )+a4_3 *F(ii,4 ) &
                    +a5_3 *F(ii,5 )+a6_3 *F(ii,6 ) &
                    +a7_3 *F(ii,7 )+a8_3 *F(ii,8 ) &
                    +a9_3 *F(ii,9 )+a10_3*F(ii,10) &
                    +a11_3*F(ii,11)
          xx(ii,4) = a1_4 *F(ii,1 )+a2_4 *F(ii,2 ) &
                    +a3_4 *F(ii,3 )+a4_4 *F(ii,4 ) &
                    +a5_4 *F(ii,5 )+a6_4 *F(ii,6 ) &
                    +a7_4 *F(ii,7 )+a8_4 *F(ii,8 ) &
                    +a9_4 *F(ii,9 )+a10_4*F(ii,10) &
                    +a11_4*F(ii,11)
          xx(ii,5) = a1_5 *F(ii,1 )+a2_5 *F(ii,2 ) &
                    +a3_5 *F(ii,3 )+a4_5 *F(ii,4 ) &
                    +a5_5 *F(ii,5 )+a6_5 *F(ii,6 ) &
                    +a7_5 *F(ii,7 )+a8_5 *F(ii,8 ) &
                    +a9_5 *F(ii,9 )+a10_5*F(ii,10) &
                    +a11_5*F(ii,11)

          xx(ii,jjmax-4) = a1_n4 *F(ii,jjmax   )+a2_n4 *F(ii,jjmax-1) &
                          +a3_n4 *F(ii,jjmax-2 )+a4_n4 *F(ii,jjmax-3) &
                          +a5_n4 *F(ii,jjmax-4 )+a6_n4 *F(ii,jjmax-5) &
                          +a7_n4 *F(ii,jjmax-6 )+a8_n4 *F(ii,jjmax-7) &
                          +a9_n4 *F(ii,jjmax-8 )+a10_n4*F(ii,jjmax-9) &
                          +a11_n4*F(ii,jjmax-10)
          xx(ii,jjmax-3) = a1_n3 *F(ii,jjmax   )+a2_n3 *F(ii,jjmax-1) &
                          +a3_n3 *F(ii,jjmax-2 )+a4_n3 *F(ii,jjmax-3) &
                          +a5_n3 *F(ii,jjmax-4 )+a6_n3 *F(ii,jjmax-5) &
                          +a7_n3 *F(ii,jjmax-6 )+a8_n3 *F(ii,jjmax-7) &
                          +a9_n3 *F(ii,jjmax-8 )+a10_n3*F(ii,jjmax-9) &
                          +a11_n3*F(ii,jjmax-10)
          xx(ii,jjmax-2) = a1_n2 *F(ii,jjmax   )+a2_n2 *F(ii,jjmax-1) &
                          +a3_n2 *F(ii,jjmax-2 )+a4_n2 *F(ii,jjmax-3) &
                          +a5_n2 *F(ii,jjmax-4 )+a6_n2 *F(ii,jjmax-5) &
                          +a7_n2 *F(ii,jjmax-6 )+a8_n2 *F(ii,jjmax-7) &
                          +a9_n2 *F(ii,jjmax-8 )+a10_n2*F(ii,jjmax-9) &
                          +a11_n2*F(ii,jjmax-10)
          xx(ii,jjmax-1) = a1_n1 *F(ii,jjmax   )+a2_n1 *F(ii,jjmax-1) &
                          +a3_n1 *F(ii,jjmax-2 )+a4_n1 *F(ii,jjmax-3) &
                          +a5_n1 *F(ii,jjmax-4 )+a6_n1 *F(ii,jjmax-5) &
                          +a7_n1 *F(ii,jjmax-6 )+a8_n1 *F(ii,jjmax-7) &
                          +a9_n1 *F(ii,jjmax-8 )+a10_n1*F(ii,jjmax-9) &
                          +a11_n1*F(ii,jjmax-10)
          xx(ii,jjmax  ) = a1_n *F(ii,jjmax   )+a2_n *F(ii,jjmax-1) &
                          +a3_n *F(ii,jjmax-2 )+a4_n *F(ii,jjmax-3) &
                          +a5_n *F(ii,jjmax-4 )+a6_n *F(ii,jjmax-5) &
                          +a7_n *F(ii,jjmax-6 )+a8_n *F(ii,jjmax-7) &
                          +a9_n *F(ii,jjmax-8 )+a10_n*F(ii,jjmax-9) &
                          +a11_n*F(ii,jjmax-10)
        enddo

      endif

    end subroutine set_rhs
!-------------------------------------------------------------------------------
    subroutine set_coefficients

      !Point 1
      !0th order, no filter
      af_1 = 0.0d0
      a1_1 = 1.0d0
      a2_1 = 0.0d0
      a3_1 = 0.0d0
      a4_1 = 0.0d0
      a5_1 = 0.0d0
      a6_1 = 0.0d0
      a7_1 = 0.0d0
      a8_1 = 0.0d0
      a9_1 = 0.0d0
      a10_1= 0.0d0
      a11_1= 0.0d0

      !Point 2
      af_2 = afs(2)
      if (icompact/=1) then
        af_2 = 0.0d0
      endif
      if (isc(2)==2) then
        !2nd order, central
        a1_2 = 1.0d0/4.0d0+af_2/2.0d0
        a2_2 = 1.0d0/2.0d0+af_2
        a3_2 = 1.0d0/4.0d0+af_2/2.0d0
        a4_2 = 0.0d0
        a5_2 = 0.0d0
        a6_2 = 0.0d0
        a7_2 = 0.0d0
        a8_2 = 0.0d0
        a9_2 = 0.0d0
        a10_2= 0.0d0
        a11_2= 0.0d0
      else if (isc(2)==4) then
        !4th order, asymmetric
        a1_2 = 1.0d0/16.0d0+7.0d0*af_2/8.0d0
        a2_2 = 3.0d0/4.0d0+af_2/2.0d0
        a3_2 = 3.0d0/8.0d0+af_2/4.0d0
        a4_2 =-1.0d0/4.0d0+af_2/2.0d0
        a5_2 = 1.0d0/16.0d0-af_2/8.0d0
        a6_2 = 0.0d0
        a7_2 = 0.0d0
        a8_2 = 0.0d0
        a9_2 = 0.0d0
        a10_2= 0.0d0
        a11_2= 0.0d0
      else if (isc(2)==6) then
        !6th order, asymmetric
        a1_2 = 1.0d0/64.0d0+31.0d0*af_2/32.0d0
        a2_2 = 29.0d0/32.0d0+3.0d0*af_2/16.0d0
        a3_2 = 15.0d0/64.0d0+17.0d0*af_2/32.0d0
        a4_2 =-5.0d0/16.0d0+5.0d0*af_2/8.0d0
        a5_2 = 15.0d0/64.0d0-15.0d0*af_2/32.0d0
        a6_2 =-3.0d0/32.0d0+3.0d0*af_2/16.0d0
        a7_2 = 1.0d0/64.0d0-af_2/32.0d0
        a8_2 = 0.0d0
        a9_2 = 0.0d0
        a10_2= 0.0d0
        a11_2= 0.0d0
      else if (isc(2)==8) then
        !8th order, asymmetric
        a1_2 = 1.0d0/256.0d0+127.0d0*af_2/128.0d0
        a2_2 = 31.0d0/32.0d0+af_2/16.0d0
        a3_2 = 7.0d0/64.0d0+25.0d0*af_2/32.0d0
        a4_2 =-7.0d0/32.0d0+7.0d0*af_2/16.0d0
        a5_2 = 7.0d0*(5.0d0-10.0d0*af_2)/128.0d0
        a6_2 =-7.0d0/32.0d0+7.0d0*af_2/16.0d0
        a7_2 = 7.0d0/64.0d0-7.0d0*af_2/32.0d0
        a8_2 =-1.0d0/32.0d0+af_2/16.0d0
        a9_2 = 1.0d0/256.0d0-af_2/128.0d0
        a10_2= 0.0d0
        a11_2= 0.0d0
      else if (isc(2)==10) then
        !10th order, asymmetric
        a1_2 = (1.0d0+1022.0d0*af_2)/1024.0d0
        a2_2 = (507.0d0+10.0d0*af_2)/512.0d0
        a3_2 = (45.0d0+934.0d0*af_2)/1024.0d0
        a4_2 = 15.0d0*(-1.0d0+2.0d0*af_2)/128.0d0
        a5_2 = 105.0d0*(1.0d0-2.0d0*af_2)/512.0d0
        a6_2 = 63.0d0*(-1.0d0+2.0d0*af_2)/256.0d0
        a7_2 = 105.0d0*(1.0d0-2.0d0*af_2)/512.0d0
        a8_2 = 15.0d0*(-1.0d0+2.0d0*af_2)/128.0d0
        a9_2 = 45.0d0*(1.0d0-2.0d0*af_2)/1024.0d0
        a10_2= 5.0d0*(-1.0d0+2.0d0*af_2)/512.0d0
        a11_2= (1.0d0-2.0d0*af_2)/1024.0d0
      else
        !0th order, no filter
        af_2 = 0.0d0
        a1_2 = 0.0d0
        a2_2 = 1.0d0
        a3_2 = 0.0d0
        a4_2 = 0.0d0
        a5_2 = 0.0d0
        a6_2 = 0.0d0
        a7_2 = 0.0d0
        a8_2 = 0.0d0
        a9_2 = 0.0d0
        a10_2= 0.0d0
        a11_2= 0.0d0
      endif

      !Point 3
      af_3 = afs(3)
      if (icompact/=1) then
        af_3 = 0.0d0
      endif
      if (isc(3)==4) then
      !4th order, central
        a1_3 =-1.0d0/16.0d0+af_3/8.0d0
        a2_3 = 1.0d0/4.0d0+af_3/2.0d0
        a3_3 = 5.0d0/8.0d0+3.0d0*af_3/4.0d0
        a4_3 = 1.0d0/4.0d0+af_3/2.0d0
        a5_3 =-1.0d0/16.0d0+af_3/8.0d0
        a6_3 = 0.0d0
        a7_3 = 0.0d0
        a8_3 = 0.0d0
        a9_3 = 0.0d0
        a10_3= 0.0d0
        a11_3= 0.0d0
      else if (isc(3)==6) then
        !6th order, asymmetric
        a1_3 =-1.0d0/64.0d0+af_3/32.0d0
        a2_3 = 3.0d0/32.0d0+13.0d0*af_3/16.0d0
        a3_3 = 49.0d0/64.0d0+15.0d0*af_3/32.0d0
        a4_3 = 5.0d0/16.0d0+3.0d0*af_3/8.0d0
        a5_3 =-15.0d0/64.0d0+15.0d0*af_3/32.0d0
        a6_3 = 3.0d0/32.0d0-3.0d0*af_3/16.0d0
        a7_3 =-1.0d0/64.0d0+af_3/32.0d0
        a8_3 = 0.0d0
        a9_3 = 0.0d0
        a10_3= 0.0d0
        a11_3= 0.0d0
      else if (isc(3)==8) then
        !8th order, asymmetric
        a1_3 =-1.0d0/256.0d0+af_3/128.0d0
        a2_3 = 1.0d0/32.0d0+15.0d0*af_3/16.0d0
        a3_3 = 57.0d0/64.0d0+7.0d0*af_3/32.0d0
        a4_3 = 7.0d0/32.0d0+9.0d0*af_3/16.0d0
        a5_3 = 7.0d0*(-5.0d0+10.0d0*af_3)/128.0d0
        a6_3 = 7.0d0/32.0d0-7.0d0*af_3/16.0d0
        a7_3 =-7.0d0/64.0d0+7.0d0*af_3/32.0d0
        a8_3 = 1.0d0/32.0d0-af_3/16.0d0
        a9_3 =-1.0d0/256.0d0+af_3/128.0d0
        a10_3= 0.0d0
        a11_3= 0.0d0
      else if (isc(3)==10) then
        !10th order, asymmetric
        a1_3 = (-1.0d0+2.0d0*af_3)/1024.0d0
        a2_3 = (5.0d0+502.0d0*af_3)/512.0d0
        a3_3 = (979.0d0+90.0d0*af_3)/1024.0d0
        a4_3 = (15.0d0+98.0d0*af_3)/128.0d0
        a5_3 = 105.0d0*(-1.0d0+2.0d0*af_3)/512.0d0
        a6_3 = 63.0d0*(1.0d0-2.0d0*af_3)/256.0d0
        a7_3 = 105.0d0*(-1.0d0+2.0d0*af_3)/512.0d0
        a8_3 = 15.0d0*(1.0d0-2.0d0*af_3)/128.0d0
        a9_3 = 45.0d0*(-1.0d0+2.0d0*af_3)/1024.0d0
        a10_3= 5.0d0*(1.0d0-2.0d0*af_3)/512.0d0
        a11_3= (-1.0d0+2.0d0*af_3)/1024.0d0
      else
        !0th order, no filter
        af_3 = 0.0d0
        a1_3 = 0.0d0
        a2_3 = 0.0d0
        a3_3 = 1.0d0
        a4_3 = 0.0d0
        a5_3 = 0.0d0
        a6_3 = 0.0d0
        a7_3 = 0.0d0
        a8_3 = 0.0d0
        a9_3 = 0.0d0
        a10_3= 0.0d0
        a11_3= 0.0d0
      endif

      !Point 4
      af_4 = afs(4)
      if (icompact/=1) then
        af_4 = 0.0d0
      endif
      if (isc(4)==6) then
        !6th order, central
        a1_4 = 1.0d0/64.0d0-af_4/32.0d0
        a2_4 =-3.0d0/32.0d0+3.0d0*af_4/16.0d0
        a3_4 = 15.0d0/64.0d0+17.0d0*af_4/32.0d0
        a4_4 = 11.0d0/16.0d0+5.0d0*af_4/8.0d0
        a5_4 = 15.0d0/64.0d0+17.0d0*af_4/32.0d0
        a6_4 =-3.0d0/32.0d0+3.0d0*af_4/16.0d0
        a7_4 = 1.0d0/64.0d0-af_4/32.0d0
        a8_4 = 0.0d0
        a9_4 = 0.0d0
        a10_4= 0.0d0
        a11_4= 0.0d0
      else if (isc(4)==8) then
        !8th order, asymmetric
        a1_4 = 1.0d0/256.0d0-af_4/128.0d0
        a2_4 =-1.0d0/32.0d0+af_4/16.0d0
        a3_4 = 7.0d0/64.0d0+25.0d0*af_4/32.0d0
        a4_4 = 25.0d0/32.0d0+7.0d0*af_4/16.0d0
        a5_4 = 35.0d0/128.0d0+29.0d0*af_4/64.0d0
        a6_4 =-7.0d0/32.0d0+7.0d0*af_4/16.0d0
        a7_4 = 7.0d0/64.0d0-7.0d0*af_4/32.0d0
        a8_4 =-1.0d0/32.0d0+af_4/16.0d0
        a9_4 = 1.0d0/256.0d0-af_4/128.0d0
        a10_4= 0.0d0
        a11_4= 0.0d0
      else if (isc(4)==10) then
        !10th order, asymmetric
        a1_4 = (1.0d0-2.0d0*af_4)/1024.0d0
        a2_4 = 5.0d0*(-1.0d0+2.0d0*af_4)/512.0d0
        a3_4 = (45.0d0+934.0d0*af_4)/1024.0d0
        a4_4 = (113.0d0+30.0d0*af_4)/128.0d0
        a5_4 = (105.0d0+302.0d0*af_4)/512.0d0
        a6_4 = 63.0d0*(-1.0d0+2.0d0*af_4)/256.0d0
        a7_4 = 105.0d0*(1.0d0-2.0d0*af_4)/512.0d0
        a8_4 = 15.0d0*(-1.0d0+2.0d0*af_4)/128.0d0
        a9_4 = 45.0d0*(1.0d0-2.0d0*af_4)/1024.0d0
        a10_4= 5.0d0*(-1.0d0+2.0d0*af_4)/512.0d0
        a11_4= (1.0d0-2.0d0*af_4)/1024.0d0
      else
        !0th order, no filter
        af_4 = 0.0d0
        a1_4 = 0.0d0
        a2_4 = 0.0d0
        a3_4 = 0.0d0
        a4_4 = 1.0d0
        a5_4 = 0.0d0
        a6_4 = 0.0d0
        a7_4 = 0.0d0
        a8_4 = 0.0d0
        a9_4 = 0.0d0
        a10_4= 0.0d0
        a11_4= 0.0d0
      endif

      !Point 5
      af_5 = afs(5)
      if (icompact/=1) then
        af_5 = 0.0d0
      endif
      if (isc(5)==8) then
        !8th order, central
        a1_5 = -1.0d0/256.0d0+af_5/128.0d0
        a2_5 = 1.0d0/32.0d0-af_5/16.0d0
        a3_5 = (-7.0d0+14.0d0*af_5)/64.0d0
        a4_5 = (7.0d0+18.0d0*af_5)/32.0d0
        a5_5 = (93.0d0+70.0d0*af_5)/128.0d0
        a6_5 = (7.0d0+18.0d0*af_5)/32.0d0
        a7_5 = (-7.0d0+14.0d0*af_5)/64.0d0
        a8_5 = 1.0d0/32.0d0-af_5/16.0d0
        a9_5 = -1.0d0/256.0d0+af_5/128.0d0
        a10_5= 0.0d0
        a11_5= 0.0d0
      else if (isc(5)==10) then
        !10th order, asymmetric
        a1_5 = (-1.0d0+2.0d0*af_5)/1024.0d0
        a2_5 = 5.0d0*(1.0d0-2.0d0*af_5)/512.0d0
        a3_5 = 45.0d0*(-1.0d0+2.0d0*af_5)/1024.0d0
        a4_5 = (15.0d0+98.0d0*af_5)/128.0d0
        a5_5 = (407.0d0+210.0d0*af_5)/512.0d0
        a6_5 = (63.0d0+130.0d0*af_5)/256.0d0
        a7_5 = 105.0d0*(-1.0d0+2.0d0*af_5)/512.0d0
        a8_5 = 15.0d0*(1.0d0-2.0d0*af_5)/128.0d0
        a9_5 = 45.0d0*(-1.0d0+2.0d0*af_5)/1024.0d0
        a10_5= 5.0d0*(1.0d0-2.0d0*af_5)/512.0d0
        a11_5= (-1.0d0+2.0d0*af_5)/1024.0d0
      else
        !0th order, no filter
        af_5 = 0.0d0
        a1_5 = 0.0d0
        a2_5 = 0.0d0
        a3_5 = 0.0d0
        a4_5 = 0.0d0
        a5_5 = 1.0d0
        a6_5 = 0.0d0
        a7_5 = 0.0d0
        a8_5 = 0.0d0
        a9_5 = 0.0d0
        a10_5= 0.0d0
        a11_5= 0.0d0
      endif

      !Interior
      af = afs(6)
      if (icompact/=1) then
        af = 0.0d0
      endif
      if (isc(6)==8) then
        !8th order, central
        a0 = (93.0d0+70.0d0*af)/128.0d0
        a1 = (7.0d0+18.0d0*af)/32.0d0
        a2 = (-7.0d0+14.0d0*af)/64.0d0
        a3 = 1.0d0/32.0d0-af/16.0d0
        a4 = -1.0d0/256.0d0+af/128.0d0
        a5 = 0.0d0
      else if (isc(6)==10) then
        !10th order, central
        a0 = (193.0d0+126.0d0*af)/256.0d0
        a1 = (105.0d0+302.0d0*af)/512.0d0
        a2 = 15.0d0*(-1.0d0+2.0d0*af)/128.0d0
        a3 = 45.0d0*(1.0d0-2.0d0*af)/1024.0d0
        a4 = 5.0d0*(-1.0d0+2.0d0*af)/512.0d0
        a5 = (1.0d0-2.0d0*af)/1024.0d0
      else
        !0th-order, no filter
        af = 0.0d0
        a0 = 1.0d0
        a1 = 0.0d0
        a2 = 0.0d0
        a3 = 0.0d0
        a4 = 0.0d0
        a5 = 0.0d0
      endif

      !Point N-4
      af_n4 = afs(7)
      if (icompact/=1) then
        af_n4 = 0.0d0
      endif
      if (isc(7)==8) then
        !8th order, central
        a1_n4 = -1.0d0/256.0d0+af_n4/128.0d0
        a2_n4 = 1.0d0/32.0d0-af_n4/16.0d0
        a3_n4 = (-7.0d0+14.0d0*af_n4)/64.0d0
        a4_n4 = (7.0d0+18.0d0*af_n4)/32.0d0
        a5_n4 = (93.0d0+70.0d0*af_n4)/128.0d0
        a6_n4 = (7.0d0+18.0d0*af_n4)/32.0d0
        a7_n4 = (-7.0d0+14.0d0*af_n4)/64.0d0
        a8_n4 = 1.0d0/32.0d0-af_n4/16.0d0
        a9_n4 = -1.0d0/256.0d0+af_n4/128.0d0
        a10_n4= 0.0d0
        a11_n4= 0.0d0
      else if (isc(7)==10) then
        !10th order, asymmetric
        a1_n4 = (-1.0d0+2.0d0*af_n4)/1024.0d0
        a2_n4 = 5.0d0*(1.0d0-2.0d0*af_n4)/512.0d0
        a3_n4 = 45.0d0*(-1.0d0+2.0d0*af_n4)/1024.0d0
        a4_n4 = (15.0d0+98.0d0*af_n4)/128.0d0
        a5_n4 = (407.0d0+210.0d0*af_n4)/512.0d0
        a6_n4 = (63.0d0+130.0d0*af_n4)/256.0d0
        a7_n4 = 105.0d0*(-1.0d0+2.0d0*af_n4)/512.0d0
        a8_n4 = 15.0d0*(1.0d0-2.0d0*af_n4)/128.0d0
        a9_n4 = 45.0d0*(-1.0d0+2.0d0*af_n4)/1024.0d0
        a10_n4= 5.0d0*(1.0d0-2.0d0*af_n4)/512.0d0
        a11_n4= (-1.0d0+2.0d0*af_n4)/1024.0d0
      else
        !0th order, no filter
        af_n4 = 0.0d0
        a1_n4 = 0.0d0
        a2_n4 = 0.0d0
        a3_n4 = 0.0d0
        a4_n4 = 0.0d0
        a5_n4 = 1.0d0
        a6_n4 = 0.0d0
        a7_n4 = 0.0d0
        a8_n4 = 0.0d0
        a9_n4 = 0.0d0
        a10_n4= 0.0d0
        a11_n4= 0.0d0
      endif

      !Point N-3
      af_n3 = afs(8)
      if (icompact/=1) then
        af_n3 = 0.0d0
      endif
      if (isc(8)==6) then
        !6th order, central
        a1_n3 = 1.0d0/64.0d0-af_n3/32.0d0
        a2_n3 =-3.0d0/32.0d0+3.0d0*af_n3/16.0d0
        a3_n3 = 15.0d0/64.0d0+17.0d0*af_n3/32.0d0
        a4_n3 = 11.0d0/16.0d0+5.0d0*af_n3/8.0d0
        a5_n3 = 15.0d0/64.0d0+17.0d0*af_n3/32.0d0
        a6_n3 =-3.0d0/32.0d0+3.0d0*af_n3/16.0d0
        a7_n3 = 1.0d0/64.0d0-af_n3/32.0d0
        a8_n3 = 0.0d0
        a9_n3 = 0.0d0
        a10_n3= 0.0d0
        a11_n3= 0.0d0
      else if (isc(8)==8) then
        !8th order, asymmetric
        a1_n3 = 1.0d0/256.0d0-af_n3/128.0d0
        a2_n3 =-1.0d0/32.0d0+af_n3/16.0d0
        a3_n3 = 7.0d0/64.0d0+25.0d0*af_n3/32.0d0
        a4_n3 = 25.0d0/32.0d0+7.0d0*af_n3/16.0d0
        a5_n3 = 35.0d0/128.0d0+29.0d0*af_n3/64.0d0
        a6_n3 =-7.0d0/32.0d0+7.0d0*af_n3/16.0d0
        a7_n3 = 7.0d0/64.0d0-7.0d0*af_n3/32.0d0
        a8_n3 =-1.0d0/32.0d0+af_n3/16.0d0
        a9_n3 = 1.0d0/256.0d0-af_n3/128.0d0
        a10_n3= 0.0d0
        a11_n3= 0.0d0
      else if (isc(8)==10) then
        !10th order, asymmetric
        a1_n3 = (1.0d0-2.0d0*af_n3)/1024.0d0
        a2_n3 = 5.0d0*(-1.0d0+2.0d0*af_n3)/512.0d0
        a3_n3 = (45.0d0+934.0d0*af_n3)/1024.0d0
        a4_n3 = (113.0d0+30.0d0*af_n3)/128.0d0
        a5_n3 = (105.0d0+302.0d0*af_n3)/512.0d0
        a6_n3 = 63.0d0*(-1.0d0+2.0d0*af_n3)/256.0d0
        a7_n3 = 105.0d0*(1.0d0-2.0d0*af_n3)/512.0d0
        a8_n3 = 15.0d0*(-1.0d0+2.0d0*af_n3)/128.0d0
        a9_n3 = 45.0d0*(1.0d0-2.0d0*af_n3)/1024.0d0
        a10_n3= 5.0d0*(-1.0d0+2.0d0*af_n3)/512.0d0
        a11_n3= (1.0d0-2.0d0*af_n3)/1024.0d0
      else
        !0th order, no filter
        af_n3 = 0.0d0
        a1_n3 = 0.0d0
        a2_n3 = 0.0d0
        a3_n3 = 0.0d0
        a4_n3 = 1.0d0
        a5_n3 = 0.0d0
        a6_n3 = 0.0d0
        a7_n3 = 0.0d0
        a8_n3 = 0.0d0
        a9_n3 = 0.0d0
        a10_n3= 0.0d0
        a11_n3= 0.0d0
      endif

      !Point N-2
      af_n2 = afs(9)
      if (icompact/=1) then
        af_n2 = 0.0d0
      endif
      if (isc(9)==4) then
        !4th order, central
        a1_n2 =-1.0d0/16.0d0+af_n2/8.0d0
        a2_n2 = 1.0d0/4.0d0+af_n2/2.0d0
        a3_n2 = 5.0d0/8.0d0+3.0d0*af_n2/4.0d0
        a4_n2 = 1.0d0/4.0d0+af_n2/2.0d0
        a5_n2 =-1.0d0/16.0d0+af_n2/8.0d0
        a6_n2 = 0.0d0
        a7_n2 = 0.0d0
        a8_n2 = 0.0d0
        a9_n2 = 0.0d0
        a10_n2= 0.0d0
        a11_n2= 0.0d0
      else if (isc(9)==6) then
        !6th order, asymmetric
        a1_n2 =-1.0d0/64.0d0+af_n2/32.0d0
        a2_n2 = 3.0d0/32.0d0+13.0d0*af_n2/16.0d0
        a3_n2 = 49.0d0/64.0d0+15.0d0*af_n2/32.0d0
        a4_n2 = 5.0d0/16.0d0+3.0d0*af_n2/8.0d0
        a5_n2 =-15.0d0/64.0d0+15.0d0*af_n2/32.0d0
        a6_n2 = 3.0d0/32.0d0-3.0d0*af_n2/16.0d0
        a7_n2 =-1.0d0/64.0d0+af_n2/32.0d0
        a8_n2 = 0.0d0
        a9_n2 = 0.0d0
        a10_n2= 0.0d0
        a11_n2= 0.0d0
      else if (isc(9)==8) then
        !8th order, asymmetric
        a1_n2 =-1.0d0/256.0d0+af_n2/128.0d0
        a2_n2 = 1.0d0/32.0d0+15.0d0*af_n2/16.0d0
        a3_n2 = 57.0d0/64.0d0+7.0d0*af_n2/32.0d0
        a4_n2 = 7.0d0/32.0d0+9.0d0*af_n2/16.0d0
        a5_n2 = 7.0d0*(-5.0d0+10.0d0*af_n2)/128.0d0
        a6_n2 = 7.0d0/32.0d0-7.0d0*af_n2/16.0d0
        a7_n2 =-7.0d0/64.0d0+7.0d0*af_n2/32.0d0
        a8_n2 = 1.0d0/32.0d0-af_n2/16.0d0
        a9_n2 =-1.0d0/256.0d0+af_n2/128.0d0
        a10_n2= 0.0d0
        a11_n2= 0.0d0
      else if (isc(9)==10) then
        !10th order, asymmetric
        a1_n2 = (-1.0d0+2.0d0*af_n2)/1024.0d0
        a2_n2 = (5.0d0+502.0d0*af_n2)/512.0d0
        a3_n2 = (979.0d0+90.0d0*af_n2)/1024.0d0
        a4_n2 = (15.0d0+98.0d0*af_n2)/128.0d0
        a5_n2 = 105.0d0*(-1.0d0+2.0d0*af_n2)/512.0d0
        a6_n2 = 63.0d0*(1.0d0-2.0d0*af_n2)/256.0d0
        a7_n2 = 105.0d0*(-1.0d0+2.0d0*af_n2)/512.0d0
        a8_n2 = 15.0d0*(1.0d0-2.0d0*af_n2)/128.0d0
        a9_n2 = 45.0d0*(-1.0d0+2.0d0*af_n2)/1024.0d0
        a10_n2= 5.0d0*(1.0d0-2.0d0*af_n2)/512.0d0
        a11_n2= (-1.0d0+2.0d0*af_n2)/1024.0d0
      else
        !0th order, no filter
        af_n2 = 0.0d0
        a1_n2 = 0.0d0
        a2_n2 = 0.0d0
        a3_n2 = 1.0d0
        a4_n2 = 0.0d0
        a5_n2 = 0.0d0
        a6_n2 = 0.0d0
        a7_n2 = 0.0d0
        a8_n2 = 0.0d0
        a9_n2 = 0.0d0
        a10_n2= 0.0d0
        a11_n2= 0.0d0
      endif

      !Point N-1
      af_n1 = afs(10)
      if (icompact/=1) then
        af_n1 = 0.0d0
      endif
      if (isc(10)==2) then
        !2nd order, central
        a1_n1 = 1.0d0/4.0d0+af_n1/2.0d0
        a2_n1 = 1.0d0/2.0d0+af_n1
        a3_n1 = 1.0d0/4.0d0+af_n1/2.0d0
        a4_n1 = 0.0d0
        a5_n1 = 0.0d0
        a6_n1 = 0.0d0
        a7_n1 = 0.0d0
        a8_n1 = 0.0d0
        a9_n1 = 0.0d0
        a10_n1= 0.0d0
        a11_n1= 0.0d0
      else if (isc(10)==4) then
        !4th order, asymmetric
        a1_n1 = 1.0d0/16.0d0+7.0d0*af_n1/8.0d0
        a2_n1 = 3.0d0/4.0d0+af_n1/2.0d0
        a3_n1 = 3.0d0/8.0d0+af_n1/4.0d0
        a4_n1 =-1.0d0/4.0d0+af_n1/2.0d0
        a5_n1 = 1.0d0/16.0d0-af_n1/8.0d0
        a6_n1 = 0.0d0
        a7_n1 = 0.0d0
        a8_n1 = 0.0d0
        a9_n1 = 0.0d0
        a10_n1= 0.0d0
        a11_n1= 0.0d0
      else if (isc(10)==6) then
        !6th order, asymmetric
        a1_n1 = 1.0d0/64.0d0+31.0d0*af_n1/32.0d0
        a2_n1 = 29.0d0/32.0d0+3.0d0*af_n1/16.0d0
        a3_n1 = 15.0d0/64.0d0+17.0d0*af_n1/32.0d0
        a4_n1 =-5.0d0/16.0d0+5.0d0*af_n1/8.0d0
        a5_n1 = 15.0d0/64.0d0-15.0d0*af_n1/32.0d0
        a6_n1 =-3.0d0/32.0d0+3.0d0*af_n1/16.0d0
        a7_n1 = 1.0d0/64.0d0-af_n1/32.0d0
        a8_n1 = 0.0d0
        a9_n1 = 0.0d0
        a10_n1= 0.0d0
        a11_n1= 0.0d0
      else if (isc(10)==8) then
        !8th order, asymmetric
        a1_n1 = 1.0d0/256.0d0+127.0d0*af_n1/128.0d0
        a2_n1 = 31.0d0/32.0d0+af_n1/16.0d0
        a3_n1 = 7.0d0/64.0d0+25.0d0*af_n1/32.0d0
        a4_n1 =-7.0d0/32.0d0+7.0d0*af_n1/16.0d0
        a5_n1 = 7.0d0*(5.0d0-10.0d0*af_n1)/128.0d0
        a6_n1 =-7.0d0/32.0d0+7.0d0*af_n1/16.0d0
        a7_n1 = 7.0d0/64.0d0-7.0d0*af_n1/32.0d0
        a8_n1 =-1.0d0/32.0d0+af_n1/16.0d0
        a9_n1 = 1.0d0/256.0d0-af_n1/128.0d0
        a10_n1= 0.0d0
        a11_n1= 0.0d0
      else if (isc(10)==10) then
        !10th order, asymmetric
        a1_n1 = (1.0d0+1022.0d0*af_n1)/1024.0d0
        a2_n1 = (507.0d0+10.0d0*af_n1)/512.0d0
        a3_n1 = (45.0d0+934.0d0*af_n1)/1024.0d0
        a4_n1 = 15.0d0*(-1.0d0+2.0d0*af_n1)/128.0d0
        a5_n1 = 105.0d0*(1.0d0-2.0d0*af_n1)/512.0d0
        a6_n1 = 63.0d0*(-1.0d0+2.0d0*af_n1)/256.0d0
        a7_n1 = 105.0d0*(1.0d0-2.0d0*af_n1)/512.0d0
        a8_n1 = 15.0d0*(-1.0d0+2.0d0*af_n1)/128.0d0
        a9_n1 = 45.0d0*(1.0d0-2.0d0*af_n1)/1024.0d0
        a10_n1= 5.0d0*(-1.0d0+2.0d0*af_n1)/512.0d0
        a11_n1= (1.0d0-2.0d0*af_n1)/1024.0d0
      else
        !0th order, no filter
        af_n1 = 0.0d0
        a1_n1 = 0.0d0
        a2_n1 = 1.0d0
        a3_n1 = 0.0d0
        a4_n1 = 0.0d0
        a5_n1 = 0.0d0
        a6_n1 = 0.0d0
        a7_n1 = 0.0d0
        a8_n1 = 0.0d0
        a9_n1 = 0.0d0
        a10_n1= 0.0d0
        a11_n1= 0.0d0
      endif

      !Point N
      !0th order, no filter
      af_n = 0.0d0
      a1_n = 1.0d0
      a2_n = 0.0d0
      a3_n = 0.0d0
      a4_n = 0.0d0
      a5_n = 0.0d0
      a6_n = 0.0d0
      a7_n = 0.0d0
      a8_n = 0.0d0
      a9_n = 0.0d0
      a10_n= 0.0d0
      a11_n= 0.0d0

    end subroutine set_coefficients

  end subroutine calc_lp_filter_align

end module mdl_lp_filter
