module mdl_input_manager
#ifdef _USE_MPI
  use mpi
#endif
  implicit none

  integer, private, parameter ::    line_length = 300
  integer, private, parameter ::     key_length = 100
  integer, private, parameter ::   value_length = 100
  integer, private, parameter :: comment_length = 100

  !Key list
  character(key_length), private, parameter :: key_1_1 = "key_decompo_fname"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_2_1 = "key_exchange_domain_halo_fname"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_3_1 = "key_steps_end_step"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_4_1 = "key_readgrid_fname"
  character(key_length), private, parameter :: key_4_2 = "key_readgrid_scale"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_5_1 = "key_restart_step"
  character(key_length), private, parameter :: key_5_2 = "key_restart_mode"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_6_1 = "key_time_delta_mode"
  character(key_length), private, parameter :: key_6_2 = "key_time_delta_cfl"
  character(key_length), private, parameter :: key_6_3 = "key_time_delta_min_max_limits"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_7_1 = "key_time_integral_mode"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_8_1 = "key_explicit_mode"
  !-----------------------------------------------------------------------------
  character(key_length), private, parameter :: key_9_1 = "key_implicit_mode"
  character(key_length), private, parameter :: key_9_2 = "key_implicit_inner_iteration"
  character(key_length), private, parameter :: key_9_3 = "key_implicit_relaxation"
  !-----------------------------------------------------------------------------

  !Value list
  character(value_length), private, save               :: decompo_fname
  !-----------------------------------------------------------------------------
  character(value_length), private, save               :: exchange_domain_halo_fname
  !-----------------------------------------------------------------------------
  integer                , private, save               :: steps_end_step = 0
  !-----------------------------------------------------------------------------
  character(value_length), private, save               :: readgrid_fname
  real(8)                , private, save               :: readgrid_scale = 1.0d0
  !-----------------------------------------------------------------------------
  integer                , private, save               :: restart_step = 0
  integer                , private, save               :: restart_mode = 0
  !-----------------------------------------------------------------------------
  integer                , private, save               :: time_delta_mode = 1
  real(8)                , private, save               :: time_delta_cfl  = 1.0d0
  real(8)                , private, save               :: time_delta_min_limit = 0.0d0
  real(8)                , private, save               :: time_delta_max_limit = 1.0d12
  !-----------------------------------------------------------------------------
  integer                , private, save               :: time_integral_mode = 0
  !-----------------------------------------------------------------------------
  integer                , private, save               :: explicit_mode = 0
  !-----------------------------------------------------------------------------
  integer                , private, save               :: implicit_mode = 0
  integer                , private, save               :: implicit_inner_iteration = 1
  real(8)                , private, save               :: implicit_relaxation = 1.0d0
  !-----------------------------------------------------------------------------

contains

  subroutine input_manager(fname)
    character(*), intent(in) :: fname

    integer :: myrank,nprocs,iproc,ierror
    integer :: line_num
    integer :: n

    character(   line_length) :: str_line
    character(    key_length) :: str_key
    character(  value_length) :: str_value
    character(comment_length) :: str_comment
    integer :: ierror_decompose

    myrank = 0
    nprocs = 1

#ifdef _USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierror)
#endif

    do iproc=0,nprocs-1
      if (myrank==iproc) then

        open(31,file=fname,form='formatted')

        read(31,*) line_num

        do n=1,line_num

          read(31,'(a)') str_line

          call decompose_line(str_line,str_key,str_value,str_comment,ierror_decompose)

          if (ierror_decompose/=0) cycle

          if (str_key==key_1_1) then
            decompo_fname = str_value
          else if (str_key==key_2_1) then
            exchange_domain_halo_fname = str_value
          else if (str_key==key_3_1) then
            read(str_value,*) steps_end_step
          else if (str_key==key_4_1) then
            readgrid_fname = str_value
          else if (str_key==key_4_2) then
            read(str_value,*) readgrid_scale
          else if (str_key==key_5_1) then
            read(str_value,*) restart_step
          else if (str_key==key_5_2) then
            read(str_value,*) restart_mode
          else if (str_key==key_6_1) then
            read(str_value,*) time_delta_mode
          else if (str_key==key_6_2) then
            read(str_value,*) time_delta_cfl
          else if (str_key==key_6_3) then
            read(str_value,*) time_delta_min_limit,time_delta_max_limit
          else if (str_key==key_7_1) then
            read(str_value,*) time_integral_mode
          else if (str_key==key_8_1) then
            read(str_value,*) explicit_mode
          else if (str_key==key_9_1) then
            read(str_value,*) implicit_mode
          else if (str_key==key_9_2) then
            read(str_value,*) implicit_inner_iteration
          else if (str_key==key_9_3) then
            read(str_value,*) implicit_relaxation
          endif

        enddo

        close(31)

      endif

#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

    enddo

    if (myrank==0) then
      write(6,*) '/-----------------------------------------------------------/'
      write(6,*) '/                      input_manager                        /'
      write(6,*) '/-----------------------------------------------------------/'

      write(6,'(A50,A3,  A50 )') key_1_1, " = ", decompo_fname
      write(6,'(A50,A3,  A50 )') key_2_1, " = ", exchange_domain_halo_fname
      write(6,'(A50,A3,  I6  )') key_3_1, " = ", steps_end_step
      write(6,'(A50,A3,  A50 )') key_4_1, " = ", readgrid_fname
      write(6,'(A50,A3, E12.5)') key_4_2, " = ", readgrid_scale
      write(6,'(A50,A3,  I6  )') key_5_1, " = ", restart_step
      write(6,'(A50,A3,  I6  )') key_5_2, " = ", restart_mode
      write(6,'(A50,A3,  I6  )') key_6_1, " = ", time_delta_mode
      write(6,'(A50,A3, E12.5)') key_6_2, " = ", time_delta_cfl
      write(6,'(A50,A3,2E12.5)') key_6_3, " = ", time_delta_min_limit,time_delta_max_limit
      write(6,'(A50,A3,  I6  )') key_7_1, " = ", time_integral_mode
      write(6,'(A50,A3,  I6  )') key_8_1, " = ", explicit_mode
      write(6,'(A50,A3,  I6  )') key_9_1, " = ", implicit_mode
      write(6,'(A50,A3,  I6  )') key_9_2, " = ", implicit_inner_iteration
      write(6,'(A50,A3, E12.5)') key_9_3, " = ", implicit_relaxation
      call flush(6)
    endif

#ifdef _USE_MPI
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

  end subroutine input_manager
  !-----------------------------------------------------------------------------
  subroutine input_decompo(fname)
    character(*), intent(out) :: fname
    fname = decompo_fname
  end subroutine input_decompo
  !-----------------------------------------------------------------------------
  subroutine input_exchange_domain_halo(fname)
    character(*), intent(out) :: fname
    fname = exchange_domain_halo_fname
  end subroutine input_exchange_domain_halo
  !-----------------------------------------------------------------------------
  subroutine input_steps(nend)
    integer, intent(out) :: nend
    nend = steps_end_step
  end subroutine input_steps
  !-----------------------------------------------------------------------------
  subroutine input_readgrid(fname,scale)
    character(*), intent(out) :: fname
    real(8)     , intent(out) :: scale
    fname = readgrid_fname
    scale = readgrid_scale
  end subroutine input_readgrid
  !-----------------------------------------------------------------------------
  subroutine input_restart(step,mode)
    integer, intent(out) :: step
    integer, intent(out) :: mode
    step = restart_step
    mode = restart_mode
  end subroutine input_restart
  !-----------------------------------------------------------------------------
  subroutine input_time_delta(mode,cfl,dt_min,dt_max)
    integer, intent(out) :: mode
    real(8), intent(out) :: cfl
    real(8), intent(out) :: dt_min
    real(8), intent(out) :: dt_max
    mode = time_delta_mode
    cfl  = time_delta_cfl
    dt_min = time_delta_min_limit
    dt_max = time_delta_max_limit
  end subroutine input_time_delta
  !-----------------------------------------------------------------------------
  subroutine input_time_integral(mode)
    integer, intent(out) :: mode
    mode = time_integral_mode
  end subroutine input_time_integral
  !-----------------------------------------------------------------------------
  subroutine input_explicit(mode)
    integer, intent(out) :: mode
    mode = explicit_mode
  end subroutine input_explicit
  !-----------------------------------------------------------------------------
  subroutine input_implicit(mode,inner_iteration,relaxation_factor)
    integer, intent(out) :: mode
    integer, intent(out) :: inner_iteration
    real(8), intent(out) :: relaxation_factor
    mode = implicit_mode
    inner_iteration = implicit_inner_iteration
    relaxation_factor = implicit_relaxation
  end subroutine input_implicit
  !-----------------------------------------------------------------------------
  subroutine decompose_line(str_line,str_key,str_value,str_comment,ierror)
    character(*), intent(in ) :: str_line
    character(*), intent(out) :: str_key
    character(*), intent(out) :: str_value
    character(*), intent(out) :: str_comment
    integer, intent(out) :: ierror

    character(   line_length) :: str_line_tmp
    character(    key_length) :: str_key_tmp
    character(  value_length) :: str_value_tmp
    character(comment_length) :: str_comment_tmp

    integer :: i,j
    integer :: num_delimiter
    integer, dimension(2) :: index_delimiter

    str_line_tmp = str_line

    index_delimiter(:) = 0
    num_delimiter = 0
    i = 1
    do
      j = index(str_line_tmp(i:),'#')
      if (j==0) exit
      num_delimiter = num_delimiter+1
      if (num_delimiter>2) exit
      index_delimiter(num_delimiter) = i+j-1
      i = i+j
    enddo

    if (num_delimiter==0) then
      ierror = 1
      return
    endif

    str_key_tmp = str_line_tmp(1:index_delimiter(1)-1)
    str_key_tmp = adjustl(str_key_tmp)
    str_key_tmp = trim(str_key_tmp)

    if (num_delimiter==1) then
      str_value_tmp = str_line_tmp(index_delimiter(1)+1:)
    else
      str_value_tmp = str_line_tmp(index_delimiter(1)+1:index_delimiter(2)-1)
    endif
    str_value_tmp = adjustl(str_value_tmp)
    str_value_tmp = trim(str_value_tmp)

    if (num_delimiter>1) then
      str_comment_tmp = str_line_tmp(index_delimiter(2)+1:)
      str_comment_tmp = adjustl(str_comment_tmp)
      str_comment_tmp = trim(str_comment_tmp)
    endif

    str_key     = str_key_tmp
    str_value   = str_value_tmp
    str_comment = str_comment_tmp

    ierror = 0

  end subroutine decompose_line

end module mdl_input_manager
