module mdl_look_up_table
  implicit none

  type look_up_table
    integer :: iallocate = 0
    integer :: nx = 0
    integer :: ny = 0
    real(8) :: xmin = 0.0d0
    real(8) :: xmax = 0.0d0
    real(8) :: ymin = 0.0d0
    real(8) :: ymax = 0.0d0
    real(8), pointer, dimension(:,:) :: val
  end type look_up_table

contains
  !-----------------------------------------------------------------------------
  subroutine look_up_table_allocate(TBL,nx,ny,xmin,xmax,ymin,ymax)
    type(look_up_table), intent(inout) :: TBL
    integer, intent(in) :: nx,ny
    real(8), intent(in) :: xmin,xmax,ymin,ymax

    if (TBL%iallocate==1) then
      deallocate(TBL%val)
      TBL%iallocate = 0
    endif

    TBL%nx = nx
    TBL%ny = ny
    TBL%xmin = xmin
    TBL%xmax = xmax
    TBL%ymin = ymin
    TBL%ymax = ymax

    allocate(TBL%val(nx,ny))
    TBL%iallocate = 1
    TBL%val = 0.0d0
  end subroutine look_up_table_allocate
  !-----------------------------------------------------------------------------
  subroutine look_up_table_deallocate(TBL)
    type(look_up_table), intent(inout) :: TBL
    if (TBL%iallocate==1) then
      deallocate(TBL%val)
      TBL%iallocate = 0
      TBL%nx = 0
      TBL%ny = 0
      TBL%xmin = 0.0d0
      TBL%xmax = 0.0d0
      TBL%ymin = 0.0d0
      TBL%ymax = 0.0d0
    endif
  end subroutine look_up_table_deallocate
  !-----------------------------------------------------------------------------
  subroutine look_up_table_calc_func_single(TBL,xarg,xref,yarg,yref,func,fref,ivalid)
    type(look_up_table), intent(in) :: TBL
    real(8), intent(in ) :: xarg,yarg
    real(8), intent(in ) :: xref,yref
    real(8), intent(out) :: func
    real(8), intent(in ) :: fref
    integer, intent(out) :: ivalid

    real(8) :: x,y,dx,dy
    integer :: ix,iy
    real(8) :: f00,f01,f10,f11,fx0,fx1
    real(8) :: xw,yw
    real(8) :: rfref

    rfref = 1.0d0/fref

    ivalid = 1

    x = xarg*xref
    y = yarg*yref

    if (x<TBL%xmin.or.TBL%xmax<x.or.y<TBL%ymin.or.TBL%ymax<y) ivalid = -1

    x = dmax1(x,TBL%xmin)
    x = dmin1(x,TBL%xmax)

    y = dmax1(y,TBL%ymin)
    y = dmin1(y,TBL%ymax)

    dx = (TBL%xmax-TBL%xmin)/dble(TBL%nx-1)
    dy = (TBL%ymax-TBL%ymin)/dble(TBL%ny-1)

    ix = int((x-TBL%xmin)/dx)+1
    xw = (x-(dx*dble(ix-1)+TBL%xmin))/dx

    iy = int((y-TBL%ymin)/dy)+1
    yw = (y-(dy*dble(iy-1)+TBL%ymin))/dy

    f00 = TBL%val(ix  ,iy  )
    f10 = TBL%val(ix+1,iy  )
    f01 = TBL%val(ix  ,iy+1)
    f11 = TBL%val(ix+1,iy+1)

    fx0 = (1.0d0-xw)*f00+xw*f10
    fx1 = (1.0d0-xw)*f01+xw*f11

    func = ((1.0d0-yw)*fx0+yw*fx1)*rfref

  end subroutine look_up_table_calc_func_single
  !-----------------------------------------------------------------------------
  subroutine look_up_table_calc_func(TBL,XARG,xref,YARG,yref,FUNC,fref, &
                                     ista,iend,jsta,jend,ksta,kend,ivalid)
    type(look_up_table), intent(in) :: TBL
    real(8), intent(in ), dimension(:,:,:) :: XARG
    real(8), intent(in ), dimension(:,:,:) :: YARG
    real(8), intent(in ) :: xref,yref
    real(8), intent(out), dimension(:,:,:) :: FUNC
    real(8), intent(in ) :: fref
    integer, intent(in ) :: ista,iend,jsta,jend,ksta,kend
    integer, intent(out) :: ivalid

    integer :: i,j,k,ii
    integer :: imax,jmax,kmax
    integer :: iimax

    real(8), allocatable, dimension(:) :: XARG_ALN,YARG_ALN,FUNC_ALN

    imax = iend-ista+1
    jmax = jend-jsta+1
    kmax = kend-ksta+1

    iimax = imax*jmax*kmax

    allocate(XARG_ALN(iimax))
    allocate(YARG_ALN(iimax))
    allocate(FUNC_ALN(iimax))

    do ii=1,iimax
      i = ista-1 + mod(ii-1,imax)+1
      j = jsta-1 + mod(ii-1,imax*jmax)/imax+1
      k = ksta-1 + (ii-1)/(imax*jmax)+1

      XARG_ALN(ii) = XARG(i,j,k)
      YARG_ALN(ii) = YARG(i,j,k)
    enddo

    call look_up_table_calc_func_align(TBL,XARG_ALN,xref,YARG_ALN,yref, &
                                       FUNC_ALN,fref,iimax,ivalid)

    do ii=1,iimax
      i = ista-1 + mod(ii-1,imax)+1
      j = jsta-1 + mod(ii-1,imax*jmax)/imax+1
      k = ksta-1 + (ii-1)/(imax*jmax)+1

      FUNC(i,j,k) = FUNC_ALN(ii)
    enddo

    deallocate(XARG_ALN)
    deallocate(YARG_ALN)
    deallocate(FUNC_ALN)
  end subroutine look_up_table_calc_func
  !-----------------------------------------------------------------------------
  subroutine look_up_table_calc_func_align(TBL,XARG_ALN,xref,YARG_ALN,yref,FUNC_ALN,fref, &
                                           iimax,ivalid)
    type(look_up_table), intent(in) :: TBL
    real(8), intent(in ), dimension(:) :: XARG_ALN
    real(8), intent(in ), dimension(:) :: YARG_ALN
    real(8), intent(in ) :: xref,yref
    real(8), intent(out), dimension(:) :: FUNC_ALN
    real(8), intent(in ) :: fref
    integer, intent(in ) :: iimax
    integer, intent(out) :: ivalid

    integer :: ii
    real(8) :: x,y,dx,dy
    integer :: ix,iy
    real(8) :: f00,f01,f10,f11,fx0,fx1
    real(8) :: xw,yw
    real(8) :: rfref

    rfref = 1.0d0/fref

    ivalid = 1

    do ii=1,iimax
      x = XARG_ALN(ii)*xref
      y = YARG_ALN(ii)*yref

      if (x<TBL%xmin.or.TBL%xmax<x.or.y<TBL%ymin.or.TBL%ymax<y) ivalid = -1

      x = dmax1(x,TBL%xmin)
      x = dmin1(x,TBL%xmax)

      y = dmax1(y,TBL%ymin)
      y = dmin1(y,TBL%ymax)

      dx = (TBL%xmax-TBL%xmin)/dble(TBL%nx-1)
      dy = (TBL%ymax-TBL%ymin)/dble(TBL%ny-1)

      ix = int((x-TBL%xmin)/dx)+1
      xw = (x-(dx*dble(ix-1)+TBL%xmin))/dx

      iy = int((y-TBL%ymin)/dy)+1
      yw = (y-(dy*dble(iy-1)+TBL%ymin))/dy

      f00 = TBL%val(ix  ,iy  )
      f10 = TBL%val(ix+1,iy  )
      f01 = TBL%val(ix  ,iy+1)
      f11 = TBL%val(ix+1,iy+1)

      fx0 = (1.0d0-xw)*f00+xw*f10
      fx1 = (1.0d0-xw)*f01+xw*f11

      FUNC_ALN(ii) = ((1.0d0-yw)*fx0+yw*fx1)*rfref
    enddo

  end subroutine look_up_table_calc_func_align
  !-----------------------------------------------------------------------------
  subroutine look_up_table_calc_func_2d(TBL,XARG,xref,YARG,yref,FUNC,fref, &
                                        ista,iend,jsta,jend,ivalid)
    type(look_up_table), intent(in) :: TBL
    real(8), intent(in ), dimension(:,:) :: XARG
    real(8), intent(in ), dimension(:,:) :: YARG
    real(8), intent(in ) :: xref,yref
    real(8), intent(out), dimension(:,:) :: FUNC
    real(8), intent(in ) :: fref
    integer, intent(in ) :: ista,iend,jsta,jend
    integer, intent(out) :: ivalid

    integer :: i,j
    real(8) :: x,y,dx,dy
    integer :: ix,iy
    real(8) :: f00,f01,f10,f11,fx0,fx1
    real(8) :: xw,yw
    real(8) :: rfref

    rfref = 1.0d0/fref

    ivalid = 1

    do j=jsta,jend
      do i=ista,iend
        x = XARG(i,j)*xref
        y = YARG(i,j)*yref

        if (x<TBL%xmin.or.TBL%xmax<x.or.y<TBL%ymin.or.TBL%ymax<y) ivalid = -1

        x = dmax1(x,TBL%xmin)
        x = dmin1(x,TBL%xmax)

        y = dmax1(y,TBL%ymin)
        y = dmin1(y,TBL%ymax)

        dx = (TBL%xmax-TBL%xmin)/dble(TBL%nx-1)
        dy = (TBL%ymax-TBL%ymin)/dble(TBL%ny-1)

        ix = int((x-TBL%xmin)/dx)+1
        xw = (x-(dx*dble(ix-1)+TBL%xmin))/dx

        iy = int((y-TBL%ymin)/dy)+1
        yw = (y-(dy*dble(iy-1)+TBL%ymin))/dy

        f00 = TBL%val(ix  ,iy  )
        f10 = TBL%val(ix+1,iy  )
        f01 = TBL%val(ix  ,iy+1)
        f11 = TBL%val(ix+1,iy+1)

        fx0 = (1.0d0-xw)*f00+xw*f10
        fx1 = (1.0d0-xw)*f01+xw*f11

        FUNC(i,j) = ((1.0d0-yw)*fx0+yw*fx1)*rfref
      enddo
    enddo

  end subroutine look_up_table_calc_func_2d
  !-----------------------------------------------------------------------------
  subroutine look_up_table_get_xy(TBL,i,j,x,y)
    type(look_up_table), intent(in) :: TBL
    integer, intent(in ) :: i,j
    real(8), intent(out) :: x,y

    real(8) :: dx,dy

    dx = (TBL%xmax-TBL%xmin)/dble(TBL%nx-1)
    dy = (TBL%ymax-TBL%ymin)/dble(TBL%ny-1)

    x = TBL%xmin+dx*dble(i-1)
    y = TBL%ymin+dy*dble(j-1)
  end subroutine look_up_table_get_xy
  !-----------------------------------------------------------------------------
  subroutine look_up_table_write(fname,TBL)
    character(*), intent(in) :: fname
    type(look_up_table), intent(inout) :: TBL

    integer :: i,j

    101 format(2i9)
    102 format(1e22.15,1e23.15)
    103 format(1e22.15,3e23.15)

    open(23,file=fname,form='formatted')

    write(23,101) TBL%nx,TBL%ny
    write(23,102) TBL%xmin,TBL%xmax
    write(23,102) TBL%ymin,TBL%ymax

    write(23,103) ((TBL%val(i,j),i=1,TBL%nx),j=1,TBL%ny)

    close(23)

  end subroutine look_up_table_write
  !-----------------------------------------------------------------------------
  subroutine look_up_table_read(fname,TBL)
    character(*), intent(in) :: fname
    type(look_up_table), intent(inout) :: TBL

    integer :: i,j

    101 format(2i9)
    102 format(1e22.15,1e23.15)
    103 format(1e22.15,3e23.15)

    open(23,file=fname,form='formatted')

    if (TBL%iallocate==1) then
      deallocate(TBL%val)
      TBL%iallocate = 0
    endif

    read(23,101) TBL%nx,TBL%ny
    read(23,102) TBL%xmin,TBL%xmax
    read(23,102) TBL%ymin,TBL%ymax

    allocate(TBL%val(TBL%nx,TBL%ny))

    read(23,103) ((TBL%val(i,j),i=1,TBL%nx),j=1,TBL%ny)

    close(23)

  end subroutine look_up_table_read

end module mdl_look_up_table
