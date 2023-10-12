module mdl_muscl
  implicit none

contains
!------------------------------------------------------------------------------!
!F : Input function
!FL: Left-interpolation
!FR: Right-interpolation
!
!    i-1          i          i+1
!     *     |     *     |     *
!         FR_i   F_i  FL_i
!
!idir
!   <1>:xsi <2>:eta <3>:zta
!
!iperi
!   <1>:periodic <else>:bounded
!
!isc
!   <1>:1st-ordeer <2>:2nd-order <3>:3rd-order
!------------------------------------------------------------------------------!
  subroutine muscl_align(F,FL,FR,isc,iimax,jjmax)
    real(8), intent(in ), dimension(:,:) :: F
    real(8), intent(out), dimension(:,:) :: FL,FR
    integer, intent(in ) :: isc
    integer, intent(in ) :: iimax,jjmax

    real(8) :: cka,ckb,eps
    real(8) :: slm,dqpp,dqpm
    integer :: ii,jj

    eps = 1.0d-6

    if (isc==1) then
      cka = 0.0d0
      ckb = 0.0d0
    else if (isc==2) then
      cka = 0.0d0
      ckb = 1.0d0
    else
      cka = 1.0d0/3.0d0
      ckb = 1.0d0
    endif

    do jj=2,jjmax-1
      do ii=1,iimax
        dqpm = F(ii,jj  )-F(ii,jj-1)
        dqpp = F(ii,jj+1)-F(ii,jj  )

        slm = (2.0d0*dqpp*dqpm+eps)/(dqpp*dqpp+dqpm*dqpm+eps)
        FL(ii,jj+1) = F(ii,jj) + 0.25d0*ckb*slm*((1.0d0-cka*slm)*dqpm+(1.0d0+cka*slm)*dqpp)
        FR(ii,jj  ) = F(ii,jj) - 0.25d0*ckb*slm*((1.0d0+cka*slm)*dqpm+(1.0d0-cka*slm)*dqpp)
      enddo
    enddo

    do ii=1,iimax
      FL(ii,1    ) = F(ii,1    )
      FL(ii,2    ) = F(ii,1    )
      FR(ii,1    ) = F(ii,1    )
      FL(ii,jjmax+1) = F(ii,jjmax)
      FR(ii,jjmax+1) = F(ii,jjmax)
      FR(ii,jjmax  ) = F(ii,jjmax)
    enddo

  end subroutine muscl_align

end module mdl_muscl
