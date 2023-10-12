module mdl_weno
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
  subroutine weno5_align(F,FL,FR,iimax,jjmax)
    real(8), intent(in ), dimension(:,:) :: F
    real(8), intent(out), dimension(:,:) :: FL,FR
    integer, intent(in ) :: iimax,jjmax

    integer :: ii,jj
    integer :: jm2,jm1,jzr,jp1,jp2
    real(8) :: eps,cfb1,cfb2
    real(8) :: fm2,fm1,fzr,fp1,fp2
    !------------------------
    real(8) :: fr0,fr1,fr2
    real(8) :: ar0,ar1,ar2
    real(8) :: br0,br1,br2
    real(8) :: wr0,wr1,wr2,rsr
    !-------------------------
    real(8) :: fl0,fl1,fl2
    real(8) :: al0,al1,al2
    real(8) :: bl0,bl1,bl2
    real(8) :: wl0,wl1,wl2,rsl
    !-------------------------

    real(8) :: rv6

    eps = 1.0d-6

    rv6 = 1.0d0/6.0d0

    cfb1 = 13.0d0/12.0d0
    cfb2 = 1.0d0/4.0d0

    do jj=3,jjmax-2
      jm2 = jj-2
      jm1 = jj-1
      jzr = jj
      jp1 = jj+1
      jp2 = jj+2
      do ii=1,iimax
        fm2 = F(ii,jm2)
        fm1 = F(ii,jm1)
        fzr = F(ii,jzr)
        fp1 = F(ii,jp1)
        fp2 = F(ii,jp2)

        fr0 = rv6*(-fm2+5.0d0*fm1+2.0d0*fzr)
        fr1 = rv6*(2.0d0*fm1+5.0d0*fzr-fp1)
        fr2 = rv6*(11.0d0*fzr-7.0d0*fp1+2.0d0*fp2)

        br0 = eps + cfb1*(fzr-2.0d0*fm1+fm2)**2 + cfb2*(3.0d0*fzr-4.0d0*fm1+fm2)**2
        br1 = eps + cfb1*(fp1-2.0d0*fzr+fm1)**2 + cfb2*(fp1-fm1)**2
        br2 = eps + cfb1*(fp2-2.0d0*fp1+fzr)**2 + cfb2*(fp2-4.0d0*fp1+3.0d0*fzr)**2

        ar0 = 0.3d0/(br0*br0)
        ar1 = 0.6d0/(br1*br1)
        ar2 = 0.1d0/(br2*br2)

        rsr = 1.0d0/(ar0+ar1+ar2)
        wr0 = ar0*rsr
        wr1 = ar1*rsr
        wr2 = ar2*rsr

        !-----------------------------------

        fl0 = rv6*(2.0d0*fm2-7.0d0*fm1+11.0d0*fzr)
        fl1 = rv6*(-fm1+5.0d0*fzr+2.0d0*fp1)
        fl2 = rv6*(2.0d0*fzr+5.0d0*fp1-fp2)

        bl0 = eps + cfb1*(fm2-2.0d0*fm1+fzr)**2 + cfb2*(fm2-4.0d0*fm1+3.0d0*fzr)**2
        bl1 = eps + cfb1*(fm1-2.0d0*fzr+fp1)**2 + cfb2*(fm1-fp1)**2
        bl2 = eps + cfb1*(fzr-2.0d0*fp1+fp2)**2 + cfb2*(3.0d0*fzr-4.0d0*fp1+fp2)**2

        al0 = 0.1d0/(bl0*bl0)
        al1 = 0.6d0/(bl1*bl1)
        al2 = 0.3d0/(bl2*bl2)

        rsl = 1.0d0/(al0+al1+al2)
        wl0 = al0*rsl
        wl1 = al1*rsl
        wl2 = al2*rsl

        FL(ii,jj+1) = wl0*fl0+wl1*fl1+wl2*fl2
        FR(ii,jj  ) = wr0*fr0+wr1*fr1+wr2*fr2
      enddo
    enddo

    do jj=2,jjmax-1,jjmax-3
      jm1 = jj-1
      jzr = jj
      jp1 = jj+1
      do ii=1,iimax
        fm1 = F(ii,jm1)
        fzr = F(ii,jzr)
        fp1 = F(ii,jp1)

        fr0 = 0.5d0*(fm1+fzr)
        fr1 = 0.5d0*(3.0d0*fzr-fp1)

        br0 = eps + (fm1-fzr)**2
        br1 = eps + (fzr-fp1)**2

        ar0 = 2.0d0/(3.0d0*br0*br0)
        ar1 = 1.0d0/(3.0d0*br1*br1)

        rsr = 1.0d0/(ar0+ar1)
        wr0 = ar0*rsr
        wr1 = ar1*rsr

        !-----------------------------------

        fl0 = 0.5d0*(-fm1+3.0d0*fzr)
        fl1 = 0.5d0*(fzr+fp1)

        bl0 = eps + (fzr-fm1)**2
        bl1 = eps + (fp1-fzr)**2

        al0 = 1.0d0/(3.0d0*bl0*bl0)
        al1 = 2.0d0/(3.0d0*bl1*bl1)

        rsl = 1.0d0/(al0+al1)
        wl0 = al0*rsl
        wl1 = al1*rsl

        FL(ii,jj+1) = wl0*fl0+wl1*fl1
        FR(ii,jj  ) = wr0*fr0+wr1*fr1
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

  end subroutine weno5_align

end module mdl_weno
