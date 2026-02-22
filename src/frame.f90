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

module atomlib_frame
  implicit none
  private
  public :: frame_t

  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension.
  integer, parameter :: ERRLEN = 128 ! Error string length.

  type :: frame_t
    real, allocatable :: coor(:,:) ! Cartesian coordinates in nm.
    real :: box(DIM,DIM) = 0.0 ! Simulation box size in nm.
    real :: time = 0.0 ! Current simulations time in ps.
    integer :: natoms = 0 ! Number of atoms.
  contains
    ! Allocation functions
    procedure :: allocate => frame_allocate
    procedure :: reallocate => frame_reallocate
    procedure :: deallocate => frame_deallocate
    procedure :: allocated => frame_allocated
    ! Assigment/copy functions
    procedure, private :: frame_append, frame_append_coor
    procedure, private :: frame_assign, frame_assign_coor
    generic :: operator(+) => frame_append, frame_append_coor
    generic :: assignment(=) => frame_assign, frame_assign_coor
    ! I/O Functions
    procedure, private :: frame_write_formatted
    generic :: write(formatted) => frame_write_formatted
  end type frame_t  

  ! Constructor interface
  interface frame_t
    module procedure :: frame_constructor
  end interface frame_t

contains

! Constructor for frame_t derived type. 
type(frame_t) function frame_constructor (natoms, time, coor, box) result (this)
  implicit none
  integer, intent(in) :: natoms ! Number of atoms.
  real, intent(in), optional :: time ! Simulation time in ps.
  real, intent(in), optional :: coor(DIM,natoms) ! Cartesian coordinates in nm.
  real, intent(in), optional :: box(DIM,DIM) ! Simulation box size in nm.

  call this%allocate(natoms)
  if (present(coor)) this%coor = coor
  if (present(time)) this%time = time
  if (present(box)) this%box = box

end function frame_constructor

! Allocate and initialize frame_t data for specified number of atoms.
subroutine frame_allocate (this, natoms, stat, errmsg)
  implicit none
  class(frame_t), intent(inout) :: this
  integer, intent(in) :: natoms ! Number of atoms.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  catch: block 

    if (natoms < 0) then
      message = "Invalid number of particles"
      status = 1
      exit catch
    end if

    if (allocated(this%coor)) deallocate (this%coor, STAT=status)
    allocate (this%coor(DIM, natoms), SOURCE=0., STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    
    this%natoms = natoms
    this%box = 0.0
    this%time = 0.0
  
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine frame_allocate

! eallocate frame_t data for specified number of atoms (bigger or smaller).
subroutine frame_reallocate (this, natoms, stat, errmsg)
  implicit none
  class(frame_t), intent(inout) :: this
  integer, intent(in) :: natoms ! Number of atoms.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  real, allocatable :: temp(:,:)
  character(ERRLEN) :: message
  integer :: i, status
  
  catch: block 

    if (natoms <= 0) then
      message = "Invalid number of particles"
      status = 1
      exit catch
    end if
    
    this%natoms = natoms

    if (allocated(this%coor)) then
      allocate (temp(DIM, natoms), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%coor, DIM=2))
      temp(:,:i) = this%coor(:,:i)
      call move_alloc (temp, this%coor)
    else
      allocate (this%coor(DIM, natoms), SOURCE=0., STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine frame_reallocate

! Deallocate frame_t data.
subroutine frame_deallocate (this, stat, errmsg)
  implicit none
  class(frame_t), intent(inout) :: this
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  catch: block

    status = 0

    if (allocated(this%coor)) deallocate (this%coor, STAT=status, ERRMSG=message)   
    if (status /= 0) exit catch

    this%natoms = 0
    this%box = 0.
    this%time = 0.
    
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine frame_deallocate

! Check if object is (properly) allocated.
function frame_allocated (this) result (out)
  implicit none
  class(frame_t), intent(in) :: this
  logical :: out

  out = allocated(this%coor)

end function frame_allocated

! Append frame_t data to existing frame_t data. Both objects must be allocated.
type(frame_t) function frame_append (this, other) result (result)
  implicit none
  class(frame_t), intent(in) :: this ! Frame data to append to.
  type(frame_t), intent(in) :: other ! Frame data to append.
  integer :: offset, natoms
  
  offset = size(this%coor, DIM=2)
  natoms = offset + size(other%coor, DIM=2)

  call result%allocate(natoms)

  result%time = this%time
  result%coor(:,:offset) = this%coor
  result%box = this%box

  result%coor(:,offset+1:) = other%coor

end function frame_append

type(frame_t) function frame_append_coor (this, coor) result (result)
  implicit none
  class(frame_t), intent(in) :: this ! Frame data to append to.
  real, intent(in) :: coor(:,:) ! Coordinate data to append.
  integer :: offset, natoms
  
  offset = size(this%coor, DIM=2)
  natoms = offset + size(coor, DIM=2)
 
  call result%allocate(natoms)

  result%time = this%time
  result%coor(:,:offset) = this%coor
  result%box = this%box

  result%coor(:,offset+1:) = coor

end function frame_append_coor 

! Assign value to frame_t data. Delete data if already allocated.
subroutine frame_assign (this, other)
  implicit none
  class(frame_t), intent(inout) :: this
  type(frame_t), intent(in) :: other
  integer :: natoms

  if (.not. other%allocated()) then
    error stop "Objects not allocated"
  end if

  natoms = size(other%coor, DIM=2)
  call this%allocate(natoms)

  this%time = other%time  
  this%coor = other%coor
  this%box = other%box  

end subroutine frame_assign

subroutine frame_assign_coor (this, coor)
  implicit none
  class(frame_t), intent(inout) :: this
  real, intent(in) :: coor(:,:)
  integer :: natoms

  natoms = size(coor, DIM=2)
  call this%allocate(natoms)

  ! this%time = ??? 
  this%coor = coor
  ! this%box = ???

end subroutine frame_assign_coor

! Formatted write for frame_t class
subroutine frame_write_formatted (this, unit, iotype, v_list, iostat, iomsg)
  implicit none
  class(frame_t), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype 
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  character(80) :: fmt
  integer :: i

  catch: block

    if (.not. this%allocated()) then
      iostat = 1
      iomsg = "Object not allocated"
      exit catch
    end if 

    if (iotype == "LISTDIRECTED") then
      write (unit, "(i0,4x,f8.4/)", IOSTAT=iostat, IOMSG=iomsg) this%natoms, this%time
      if (iostat /= 0) exit catch

      do i = 1, size(this%coor, DIM=2)
        write (unit, "(3(2x,f8.4),/)", IOSTAT=iostat, IOMSG=iomsg) this%coor(:,i)
        if (iostat /= 0) exit catch
      end do

      write (unit, "(*(x,f8.4),/)", IOSTAT=iostat, IOMSG=iomsg) this%box(:,:)
      if (iostat /= 0) exit catch

      write (unit, "(/)", IOSTAT=iostat, IOMSG=iomsg)
      if (iostat /= 0) exit catch
    
    else if (iotype == "DT") then
      if (size(v_list) == 0) then
        fmt = "f8.4"
      else if (size(v_list) == 2) then
        write (fmt,"(a,i0,a,i0)") "f", v_list(1), ".", v_list(2)
      else
        iostat = 1
        iomsg = "integer-list for DT descriptor must contain two integers"
        exit catch
      end if
      
      write (unit, "(i0,4x,"//trim(fmt)//",/)", IOSTAT=iostat, IOMSG=iomsg) this%natoms, this%time
      if (iostat /= 0) exit catch

      do i = 1, size(this%coor, DIM=2)
        write (unit, "(3(2x,"//trim(fmt)//"),/)", IOSTAT=iostat, IOMSG=iomsg) this%coor(:,i)
        if (iostat /= 0) exit catch
      end do 

      write (unit, "(*(x,"//trim(fmt)//"),/)", IOSTAT=iostat, IOMSG=iomsg) this%box
      if (iostat /= 0) exit catch

      write (unit, "(/)", IOSTAT=iostat, IOMSG=iomsg)
      if (iostat /= 0) exit catch

    else
      iostat = 1
      iomsg = "Unsupported iotype"
      exit catch

    end if

  end block catch

end subroutine frame_write_formatted

end module atomlib_frame
