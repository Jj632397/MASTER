module mdl_sgs_wale
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  implicit none

contains

  subroutine sgs_wale(BLK)
    type(block), intent(inout) :: BLK

    real(8), parameter :: cw = 0.575d0 !0.55 <= cw <= 0.60
    real(8), parameter :: eps= 1.0d-32

    integer :: i,j,k
    real(8) :: rh
    real(8) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(8) :: s11,s22,s33,s12,s13,s23
    real(8) :: o12,o13,o23
    real(8) :: s2,o2
    real(8) :: ss11,ss22,ss33,ss12,ss13,ss23
    real(8) :: oo11,oo22,oo33,oo12,oo13,oo23
    real(8) :: vt,vi,sd,dl,cd,ar1,ar2,ar3

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ux = BLK%DULX(i,j,k)
          uy = BLK%DULY(i,j,k)
          uz = BLK%DULZ(i,j,k)
          vx = BLK%DVLX(i,j,k)
          vy = BLK%DVLY(i,j,k)
          vz = BLK%DVLZ(i,j,k)
          wx = BLK%DWLX(i,j,k)
          wy = BLK%DWLY(i,j,k)
          wz = BLK%DWLZ(i,j,k)

          rh = BLK%DNST(i,j,k)

          dl = BLK%GRSC(i,j,k)

          s11 = ux
          s22 = vy
          s33 = wz
          s12 = 0.5d0*(uy+vx)
          s13 = 0.5d0*(uz+wx)
          s23 = 0.5d0*(vz+wy)

          o12 = 0.5d0*(uy-vx)
          o13 = 0.5d0*(uz-wx)
          o23 = 0.5d0*(vz-wy)

          s2 = s11*s11+s22*s22+s33*s33+2.0d0*(s12*s12+s13*s13+s23*s23)
          o2 = 2.0d0*(o12*o12+o13*o13+o23*o23)

          ss11 = s11*s11+s12*s12+s13*s13
          ss22 = s12*s12+s22*s22+s23*s23
          ss33 = s13*s13+s23*s23+s33*s33
          ss12 = s11*s12+s12*s22+s13*s23
          ss13 = s11*s13+s12*s23+s13*s33
          ss23 = s12*s13+s22*s23+s23*s33

          oo11 =-o12*o12-o13*o13
          oo22 =-o12*o12-o23*o23
          oo33 =-o13*o13-o23*o23
          oo12 =-o13*o23
          oo13 = o12*o23
          oo23 =-o12*o13

          vi = ss11*oo11+ss22*oo22+ss33*oo33+2.0d0*(ss12*oo12+ss13*oo13+ss23*oo23)

          sd = (s2*s2+o2*o2+4.0d0*s2*o2+12.0d0*vi)/6.0d0
          sd = dmax1(sd,eps)

          cd = cw*cw*dl*dl
          ar1= sd**1.5d0
          ar2= s2**2.5d0
          ar3= sd**1.25d0

          vt = cd*ar1/(ar2+ar3)

          BLK%VSST(i,j,k) = rh*vt*RYNLD
          BLK%CKPT(i,j,k) = BLK%VSST(i,j,k)*BLK%SHCP(i,j,k)/tprntl
        enddo
      enddo
    enddo

  end subroutine sgs_wale

end module mdl_sgs_wale
