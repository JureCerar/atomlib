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

module atomlib_ndxio
  implicit none
  private
  public :: ndx_t, group_t

  ! Global constants
  integer, parameter :: MAXLEN = 128
  integer, parameter :: ERRLEN = 128

  type group_t
    character(MAXLEN) :: title = "" ! Title of the group
    integer, allocatable :: loc(:) ! Index location
    integer :: natoms = 0 ! Size of index group
  contains
    ! Allocation functions
    procedure :: allocate => group_allocate
    procedure :: reallocate => group_reallocate
    procedure :: deallocate => group_deallocate
    procedure :: allocated => group_allocated
    ! Formated I/O
    procedure, private :: group_write_formatted
    generic :: write(formatted) => group_write_formatted
  end type group_t

  type ndx_t
    type(group_t), allocatable :: group(:) ! Index group
    integer :: ngroups = 0 ! Number of groups
  contains
    ! Allocation functions
    procedure :: allocate => ndx_allocate
    procedure :: reallocate => ndx_reallocate
    procedure :: deallocate => ndx_deallocate
    procedure :: allocated => ndx_allocated
    ! Formated I/O
    procedure, private :: ndx_write_formatted
    generic :: write(formatted) => ndx_write_formatted
    ! I/O Functions
    procedure :: save => ndx_save
    procedure :: load => ndx_load
  end type ndx_t

  interface group_t
    module procedure :: group_constructor
  end interface group_t

contains

! Util: Count number of tokens in a string seperated by delimiter.
integer function cnttok (string, del)
  implicit none
  character(*), intent(in) :: string
  character(*), intent(in) :: del
  integer :: i

  cnttok = 0
  i = len_trim(string)
  do while (i > 0)
    cnttok = cnttok + 1
    i = index(trim(string(:i)), del, BACK=.true.)
  end do

end function cnttok

! Util: Count number of atoms (indeces) in one index group.
integer function indexCount (unit, iostat, iomsg) result (cnt)
  implicit none
  integer, intent(in) :: unit
  integer, intent(out) :: iostat ! Error status code. Returns zero if no error.
  character(ERRLEN), intent(out) :: iomsg ! Error message.
  character(MAXLEN) :: buffer
  integer :: i, lines

  ! Count number of atoms in group
  cnt = 0
  lines = 0

  while: do while (.True.)
    lines = lines + 1
    read (unit, "(a)", IOSTAT=iostat, IOMSG=iomsg) buffer
    if (iostat /= 0) exit while
    ! Ops ... start of next group. Go back.
    if (index(buffer, "[") /= 0) exit while
    cnt = cnt + cnttok(adjustl(buffer), " ")

  end do while

  ! Rewind by number of lines
  do i = 1, lines
    backspace (unit, IOSTAT=iostat, IOMSG=iomsg)
    if (iostat /= 0) return
  end do

end function indexCount

! Constructor for group_t derived type.
type(group_t) function group_constructor (natoms, title, loc) result (this)
  implicit none
  integer, intent(in) :: natoms ! Number of indeces
  character(*), intent(in), optional :: title ! Title of the group
  integer, intent(in), optional :: loc(natoms) ! Index data

  call this%allocate(natoms)
  if (present(title)) this%title = title
  if (present(loc)) this%loc = loc

end function group_constructor

! Allocate and initialize group_t data for specified number of atoms.
subroutine group_allocate (this, natoms, stat, errmsg)
  implicit none
  class(group_t), intent(inout) :: this
  integer, intent(in) :: natoms ! Number of atoms
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
    
    this%natoms = natoms
  
    if (allocated(this%loc)) deallocate (this%loc, STAT=status)
    allocate (this%loc(natoms), SOURCE=0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)  

end subroutine group_allocate

! Reallocate group_t data for specified number of atoms (bigger or smaller).
subroutine group_reallocate (this, natoms, stat, errmsg)
  implicit none
  class(group_t), intent(inout) :: this
  integer, intent(in) :: natoms ! New number of atoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  integer, allocatable :: temp(:)
  character(ERRLEN) :: message
  integer :: i, status
  
  catch: block 

    if (natoms <= 0) then
      message = "Invalid number of particles"
      status = 1
      exit catch
    end if

    this%natoms = natoms
    
    if (allocated(this%loc)) then
      allocate (temp(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(size(temp), size(this%loc))
      temp(:i) = this%loc(:i)
      call move_alloc (temp, this%loc)
    else
      allocate (this%loc(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine group_reallocate

! Deallocate group_t data.
subroutine group_deallocate (this, stat, errmsg)
  implicit none
  class(group_t), intent(inout) :: this
  integer, intent(out), optional :: stat
  character(*), intent(out), optional :: errmsg
  character(ERRLEN) :: message
  integer :: status

  catch: block

    this%natoms = 0

    if (allocated(this%loc)) then
      deallocate (this%loc, STAT=status, ERRMSG=message)   
      if (status /= 0) exit catch
    else
      status = 0
    end if 

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine group_deallocate

! Check if group_t is (properly) allocated.
function group_allocated (this) result (out)
  implicit none
  logical :: out
  class(group_t), intent(in) :: this

  out = allocated(this%loc)

end function group_allocated 

! Formatted write for group_t class
subroutine group_write_formatted (this, unit, iotype, v_list, iostat, iomsg)
  implicit none
  class(group_t), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype 
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  character(80) :: fmt

  catch: block

    if (iotype == "LISTDIRECTED") then
      ! Group title
      write (unit, "(a,/)", IOSTAT=iostat, IOMSG=iomsg) "[ " // trim(this%title) // " ]"
      if (iostat /= 0) exit catch

      ! Indeces in group of 15
      if (allocated(this%loc)) then
        write (unit, "(15(i0:x))", IOSTAT=iostat, IOMSG=iomsg) this%loc
        if (iostat /= 0) exit catch
        write (unit, "(/)", IOSTAT=iostat, IOMSG=iomsg) ! New line
        if (iostat /= 0) exit catch
      end if

    else if (iotype == "DT") then
      if (size(v_list) == 0) then
        fmt = "i0"
      else if (size(v_list) == 1) then
        write (fmt,"(a,i0)") "i", v_list(1)
      else
        iomsg = "Integer-list for DT descriptor must contain one integer"
        iostat = 1
        exit catch
      end if

      ! Group title
      write (unit, "(a,/)", IOSTAT=iostat, IOMSG=iomsg) "[ " // trim(this%title) // " ]"
      if (iostat /= 0) exit catch
  
      ! Indeces in group of 15
      if (allocated(this%loc)) then
        write (unit, "(15("//trim(fmt)//":x))", IOSTAT=iostat, IOMSG=iomsg) this%loc
        if (iostat /= 0) exit catch
      end if

    else
      iostat = 1
      iomsg = "Unsupported iotype"
      exit catch

    end if

  end block catch

end subroutine group_write_formatted

! ndx_t ------------------------------------------------

! Allocate and initialize ndx_t data for specified number of groups.
subroutine ndx_allocate (this, ngroups, stat, errmsg)
  implicit none
  class(ndx_t), intent(inout) :: this
  integer, intent(in) :: ngroups ! Number of index groups
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  catch: block 

    if (ngroups < 0) then
      message = "Invalid number of groups"
      status = 1
      exit catch
    end if
    
    this%ngroups = ngroups
  
    if (allocated(this%group)) deallocate (this%group, STAT=status)
    allocate (this%group(ngroups), STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)  

end subroutine ndx_allocate

! Reallocate ndx_t data for specified number of groups (bigger or smaller).
subroutine ndx_reallocate (this, ngroups, stat, errmsg)
  implicit none
  class(ndx_t), intent(inout) :: this
  integer, intent(in) :: ngroups ! New number of groups.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  type(group_t), allocatable :: temp(:)
  character(ERRLEN) :: message
  integer :: i, status
  
  catch: block 

    if (ngroups <= 0) then
      message = "Invalid number of groups"
      status = 1
      exit catch
    end if

    this%ngroups = ngroups
    
    if (allocated(this%group)) then
      allocate (temp(ngroups), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(ngroups, size(this%group))
      temp(:i) = this%group(:i)
      call move_alloc (temp, this%group)
    else
      allocate (this%group(ngroups), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine ndx_reallocate

! Deallocate ndx_t data.
subroutine ndx_deallocate (this, stat, errmsg)
  implicit none
  class(ndx_t), intent(inout) :: this
  integer, intent(out), optional :: stat
  character(*), intent(out), optional :: errmsg
  character(ERRLEN) :: message
  integer :: status

  catch: block

    this%ngroups = 0
  
    if (allocated(this%group)) then
      deallocate (this%group, STAT=status, ERRMSG=message)   
      if (status /= 0) exit catch
    else
      status = 0
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine ndx_deallocate

! Check if ndx_t is (properly) allocated.
function ndx_allocated (this) result (out)
  implicit none
  logical :: out
  class(ndx_t), intent(in) :: this

  out = allocated(this%group)

end function ndx_allocated 

! Formatted write for ndx_t class
subroutine ndx_write_formatted (this, unit, iotype, v_list, iostat, iomsg)
  implicit none
  class(ndx_t), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype 
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  character(80) :: fmt

  catch: block

    if (iotype == "LISTDIRECTED") then
      if (allocated(this%group)) then
        write (unit, *, IOSTAT=iostat, IOMSG=iomsg) this%group
        if (iostat /= 0) error stop iomsg
      end if

    else if (iotype == "DT") then
      if (size(v_list) == 0) then
        fmt = "(DT)"
      else if (size(v_list) == 1) then
        write (fmt,"(a,i0,a,i0,a)") "(DT(", v_list(1), "))"
      else
        iomsg = "Integer-list for DT descriptor must contain one integer"
        iostat = 1
        exit catch
      end if

      if (allocated(this%group)) then
        write (unit, fmt, IOSTAT=iostat, IOMSG=iomsg) this%group
        if (iostat /= 0) error stop iomsg
      end if

    else
      iostat = 1
      iomsg = "Unsupported iotype"
      exit catch

    end if

  end block catch

end subroutine ndx_write_formatted

! Save (write) index file. 
subroutine ndx_save (this, file, stat, errmsg)
  implicit none
  class(ndx_t), intent(inout) :: this
  character(*), intent(in) :: file ! Output file.
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

    open (NEWUNIT=unit, FILE=file, ACTION="write", STATUS="unknown", IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    do i = 1, size(this%group)
      if (this%group(i)%allocated()) then
        write (unit, "(a)", IOSTAT=status, IOMSG=message) "[ " // trim(this%group(i)%title) // " ]"
        if (status /= 0) exit catch
        write (unit, "(15(i0:x))", IOSTAT=status, IOMSG=message) this%group(i)%loc
        if (status /= 0) exit catch
      end if
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

end subroutine ndx_save

! Load (read) index file.
subroutine ndx_load (this, file, stat, errmsg)
  implicit none
  class(ndx_t), intent(inout) :: this
  character(*), intent(in) :: file
  integer, intent(out), optional :: stat
  character(*), intent(out), optional :: errmsg
  character(ERRLEN) :: message
  character(MAXLEN) :: buffer, title
  integer :: i, left, right, unit, ngroups, natoms, status

  catch: block

    open (NEWUNIT=unit, FILE=file, ACTION="read", STATUS="old", IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ngroups = 0
    do while (.True.) 
      read (unit, "(a)", IOSTAT=status) buffer
      if (status /= 0) exit
      if (index(buffer, "[") /= 0) ngroups = ngroups + 1
    end do
    rewind (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch
    
    if (ngroups < 1) then
      message = "File does not contain any index groups"
      status = 1
      exit catch
    end if

    call this%allocate(ngroups, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    do i = 1, ngroups
      read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
      if (status /= 0) exit catch
      left = index(buffer, "[")
      right = index(buffer, "]")
      if (left == 0 .or. right == 0) then
        message = "Cannot read index group title"
        status = 1
        exit catch
      end if
      title = trim(adjustl(buffer(left+1:right-1)))

      natoms = indexCount(unit, status, message)
      if (status /= 0) exit catch
      
      if (natoms > 0) then 
        call this%group(i)%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        this%group(i)%title = trim(title)
        read (unit, *, IOSTAT=status, IOMSG=message) this%group(i)%loc
        if (status /= 0) exit catch
      end if      
      
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

end subroutine ndx_load

end module atomlib_ndxio
