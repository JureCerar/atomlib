! This file is part of atomlib
! https://github.com/JureCerar/atomlib
!
! Copyright (C) 2022 Jure Cerar
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module atomlib_pdhio
  implicit none
  private
  public :: pdh_t, PDH_NUM_CHAR_CONST, PDH_NUM_INT_CONST, PDH_NUM_REAL_CONST

  ! Global constants
  integer, parameter :: ERRLEN = 128
  integer, parameter :: TEXT_LEN = 80
  integer, parameter :: KEY_LEN = 4

  ! Global constants
  integer, parameter :: PDH_NUM_CHAR_CONST = 16
  integer, parameter :: PDH_NUM_INT_CONST = 8
  integer, parameter :: PDH_NUM_REAL_CONST = 10

  ! .PDH file format
  type pdh_t
    character(TEXT_LEN) :: text = "" ! Title text
    character(KEY_LEN) :: key_words(16) = "" ! Key words
    integer :: int_const(8) = 0 ! Integer constants
    real :: real_const(10) = 0.0 ! Real constants
    integer :: num_points = 0 ! Num. of points 
    real, allocatable :: x(:), y(:), y_error(:) ! x, y and error
  contains
    ! Allocation functions
    procedure :: allocate => pdh_allocate
    procedure :: reallocate => pdh_reallocate
    procedure :: deallocate => pdh_deallocate
    procedure :: allocated => pdh_allocated
    ! Formatted I/O
    procedure, private :: pdh_write_formatted
    generic :: write(formatted) => pdh_write_formatted
    ! File I/O
    procedure :: save => pdh_save
    procedure :: load => pdh_load
  end type pdh_t

  ! Constructor
  interface pdh_t
    module procedure :: pdh_constructor
  end interface pdh_t

contains

! Constructor for frame_t derived type. 
type(pdh_t) function pdh_constructor (num_points, text, key_words, int_const, real_const, x, y, y_error) result (this)
  implicit none
  integer, intent(in) :: num_points ! Number of points.
  character(*), intent(in), optional :: text ! Title text.
  character(*), intent(in), optional :: key_words(16) ! Character keywords.
  integer, intent(in), optional :: int_const(8) ! Integer constant.
  real, intent(in), optional :: real_const(10) ! Integer constant.
  real, intent(in), optional :: x(num_points), y(num_points), y_error(num_points) ! x, y, and y-error data

  call this%allocate(num_points)
  if (present(text)) this%text = text
  if (present(key_words)) this%key_words = key_words
  if (present(int_const)) this%int_const = int_const
  if (present(real_const)) this%real_const = real_const
  if (present(x)) this%x = x
  if (present(y)) this%y = y
  if (present(y_error)) this%y_error = y_error

end function pdh_constructor


! Deallocate and allocate .PDH file.
subroutine pdh_allocate (this, num_points, stat, errmsg)
  implicit none
  class(pdh_t), intent(inout) :: this
  integer, intent(in) :: num_points ! Number of points.
  integer, intent(out), optional :: stat 
  character(*), intent(out), optional :: errmsg
  character(ERRLEN) :: message
  integer :: status

  catch: block

    if (num_points < 0) then
      message = "Invalid number of points"
      status = 1
      exit catch
    end if

    this%num_points = num_points  

    if (allocated(this%x)) deallocate (this%x, STAT=status, ERRMSG=message)
    allocate (this%x(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%y)) deallocate (this%y, STAT=status, ERRMSG=message)
    allocate (this%y(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%y_error)) deallocate (this%y_error, STAT=status, ERRMSG=message)
    allocate (this%y_error(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdh_allocate

! Reallocate object to new size (smaller or bigger).
subroutine pdh_reallocate (this, num_points, stat, errmsg)
  implicit none
  class(pdh_t), intent(inout) :: this
  integer, intent(in) :: num_points ! Number of points
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  real, allocatable :: temp(:)
  character(ERRLEN) :: message
  integer :: i, status

  catch: block

    if (num_points < 0) then
      message = "Invalid number of points"
      status = 1
      exit catch
    end if

    this%num_points = num_points

    if (allocated(this%x)) then
      allocate (temp(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(num_points, size(this%x))
      temp(:i) = this%x(:i)
      call move_alloc (temp, this%x)
    else
      allocate (this%x(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%y)) then
      allocate (temp(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(num_points, size(this%y))
      temp(:i) = this%y(:i)
      call move_alloc (temp, this%y)
    else
      allocate (this%y(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%y)) then
      allocate (temp(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(num_points, size(this%y_error))
      temp(:i) = this%y_error(:i)
      call move_alloc (temp, this%y_error)
    else
      allocate (this%y_error(num_points), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdh_reallocate

! Deallocate object.
subroutine pdh_deallocate (this, stat, errmsg)
  implicit none
  class(pdh_t), intent(inout) :: this
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  catch: block
  
    this%num_points = 0

    status = 0
    if (allocated(this%x)) deallocate (this%x, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    if (allocated(this%y)) deallocate (this%y, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    if (allocated(this%y_error)) deallocate (this%y_error, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdh_deallocate

! Check if object is (properly) allocated.
function pdh_allocated (this) result (out)
  implicit none
  class(pdh_t), intent(inout) :: this
  logical :: out

  out = allocated(this%x)
  out = out .and. allocated(this%y)
  out = out .and. allocated(this%y_error)

end function pdh_allocated

! Write object data.
subroutine pdh_write_formatted (this, unit, iotype, v_list, iostat, iomsg)
  implicit none
  class(pdh_t), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype 
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  character(80) :: fmt
  integer :: i

  catch: block

    if (iotype == "LISTDIRECTED") then  
      write (unit, "(a80,/)", IOSTAT=iostat, IOMSG=iomsg) this%text
      if (iostat /= 0) exit catch
      write (unit, "(16(a4,1x),/)", IOSTAT=iostat, IOMSG=iomsg) this%key_words
      if (iostat /= 0) exit catch
      write (unit, "(8(i9,1x),/)", IOSTAT=iostat, IOMSG=iomsg) this%num_points, this%int_const(2:)
      if (iostat /= 0) exit catch
      write (unit, "(2(5(e14.6,1x),/))", IOSTAT=iostat, IOMSG=iomsg) this%real_const
      if (iostat /= 0) exit catch
      do i = 1, size(this%x)
        write (unit, "(3(1pe14.6,1x),/)", IOSTAT=iostat, IOMSG=iomsg) this%x(i), this%y(i), this%y_error(i)
        if (iostat /= 0) exit catch
      end do      

    else if (iotype == "DT") then
      if (size(v_list) == 0) then
        fmt = "14.6"
      else if (size(v_list) == 2) then
        write (fmt,"(i0,a,i0)") v_list(1), ".", v_list(2)
      else
        iostat = 1
        iomsg = "integer-list for DT descriptor must contain two integers"
        exit catch
      end if
      write (unit, "(a80,/)", IOSTAT=iostat, IOMSG=iomsg) this%text
      if (iostat /= 0) exit catch
      write (unit, "(16(a4,1x),/)", IOSTAT=iostat, IOMSG=iomsg) this%key_words
      if (iostat /= 0) exit catch
      write (unit, "(8(i9,1x),/)", IOSTAT=iostat, IOMSG=iomsg) this%num_points, this%int_const(2:)
      if (iostat /= 0) exit catch
      write (unit, "(2(5(e"//trim(fmt)//",1x),/))", IOSTAT=iostat, IOMSG=iomsg) this%real_const
      if (iostat /= 0) exit catch
      do i = 1, size(this%x)
        write (unit, "(3(1pe"//trim(fmt)//",1x),/)", IOSTAT=iostat, IOMSG=iomsg) this%x(i), this%y(i), this%y_error(i)
        if (iostat /= 0) exit catch
      end do            

    else
      iostat = 1
      iomsg = "Unsupported iotype"
      exit catch

    end if

  end block catch

end subroutine pdh_write_formatted

! write .PDH to file or unit
subroutine pdh_save (this, file, stat, errmsg)
  implicit none
  class(pdh_t), intent(inout) :: this
  character(*), intent(in) :: file ! Output file name.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: i, unit, status

  catch: block
    if (.not. this%allocated()) then
      message = "Object not allocated"
      status = 1
      exit catch
    end if

    open (NEWUNIT=unit, FILE=trim(file), STATUS="unknown", ACTION="write", IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    write (unit, "(a80)", IOSTAT=status, IOMSG=message) this%text
    if (status /= 0) exit catch

    write (unit, "(16(a4,1x))", IOSTAT=status, IOMSG=message) this%key_words
    if (status /= 0) exit catch

    write (unit, "(8(i9,1x))", IOSTAT=status, IOMSG=message) this%num_points, this%int_const(2:)
    if (status /= 0) exit catch

    write (unit, "(5(e14.6,1x))", IOSTAT=status, IOMSG=message) this%real_const
    if (status /= 0) exit catch

    do i = 1, size(this%x)
      write (unit, "(3(1pe14.6,1x))", IOSTAT=status, IOMSG=message) this%x(i), this%y(i), this%y_error(i)
      if (status /= 0) exit catch
    end do

    close (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdh_save

! Read .PDH file
subroutine pdh_load (this, file, stat, errmsg)
  implicit none
  class(pdh_t), intent(inout) :: this
  character(*), intent(in) :: file ! Input file name.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: i, unit, num_points, status

  catch: block
    open (NEWUNIT=unit, FILE=trim(file), STATUS="old", ACTION="read", IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    read (unit, "(a80)", IOSTAT=status, IOMSG=message) this%text
    if (status /= 0) exit catch

    read (unit, "(16(a4,1x))", IOSTAT=status, IOMSG=message) this%key_words
    if (status /= 0) exit catch

    ! NOTE: `int_const(1)` is defined as `num_points`.
    read (unit, "(8(i9,1x))", IOSTAT=status, IOMSG=message) num_points, this%int_const(2:)
    if (status /= 0) exit catch
    this%int_const(1) = num_points

    read (unit, "(5(e14.6,1x))", IOSTAT=status, IOMSG=message) this%real_const
    if (status /= 0) exit catch

    call this%allocate(num_points, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    do i = 1, num_points
      read (unit, "(3(1pe14.6,1x))", IOSTAT=status, IOMSG=message) this%x(i), this%y(i), this%y_error(i)
      if (status /= 0) exit catch
    end do

    close (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdh_load

end module atomlib_pdhio

! program main
!   use atomlib_pdhio
!   implicit none
!   character(*), parameter :: file = "../tests/files/data.pdh"
!   type(pdh_t) :: pdh, cpy

!   call pdh%load(file)
!   cpy = pdh / 2.


!   print *, cpy

! end program main