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

program main
  use atomlib_ndxio
  implicit none
  real, parameter :: DELTA = 0.001
  type(ndx_t) :: ndx
  character(128) :: arg, message
  integer :: i, status

  ! First check if all internal functions work properly.
  call group_test ()
  call ndx_test ()

  ! Check for each file type
  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    write (*,*) "Processing file: '" // trim(arg) //"'"

    ! Load file and check values
    call ndx%load (arg)
    if (ndx%ngroups /= 5) error stop 1
    if (ndx%group(1)%title /= "System") error stop 2
    if (ndx%group(1)%natoms /= 44) error stop 3
    if (ndx%group(1)%loc(1) /= 1) error stop 4

    ! Load file /w error handling.
    call ndx%load (arg, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message

    ! Save file
    call ndx%save ("out.ndx")

    ! Save file w/ error handling
    call ndx%save ("out.ndx", STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
  end do

contains

! Test all internal functions of group_t object.
subroutine group_test ()
  implicit none
  integer, parameter :: NP = 100
  type(group_t) :: group
  character(128) :: message
  integer :: i, status

  ! Allocate reallocate and check size
  call group%allocate (NP)
  call group%reallocate (2*NP)
  call group%reallocate (NP/2)
  call group%deallocate ()

  ! Allocate reallocate and check size
  call group%allocate (NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call group%reallocate (2*NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call group%reallocate (NP/2, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call group%deallocate (STAT=status, ERRMSG=message)
  if (status /= 0) error stop

  ! Construct a dummy configuration
  call group%allocate(NP)
  group%title = "Group"
  group%loc = [(i, i = 1, NP)]

  ! Construct a dummy configuration
  group = group_t(NP, TITLE="Group", LOC=[(i, i = 1, NP)])

  ! Try derived type write functions
  write (*, *, IOSTAT=status, IOMSG=message) group
  if (status /= 0) error stop message
  write (*,"(DT)", IOSTAT=status, IOMSG=message) group
  if (status /= 0) error stop message
  write (*,"(DT(10))", IOSTAT=status, IOMSG=message) group
  if (status /= 0) error stop message

end subroutine group_test 

! Test all internal functions of ndx_t object.
subroutine ndx_test ()
  implicit none
  integer, parameter :: NG = 4, NP = 100 
  type(ndx_t) :: ndx
  character(128) :: message, buffer
  integer :: i, status

  ! Allocate reallocate and check size
  call ndx%allocate (NG)
  call ndx%reallocate (2*NG)
  call ndx%reallocate (NG/2)
  call ndx%deallocate ()

  ! Allocate reallocate and check size
  call ndx%allocate (NG, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call ndx%reallocate (2*NG, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call ndx%reallocate (NG/2, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call ndx%deallocate (STAT=status, ERRMSG=message)
  if (status /= 0) error stop

  ! Construct a dummy configuration.
  call ndx%allocate(NG)
  do i = 1, ndx%ngroups
    write (buffer, "(a,x,i0)") "Group", i
    ndx%group(i) = group_t(NP, buffer, [(i, i = 1, NP)])
  end do

  ! Try derived type write functions
  write (*, *, IOSTAT=status, IOMSG=message) ndx
  if (status /= 0) error stop message
  write (*,"(DT)", IOSTAT=status, IOMSG=message) ndx
  if (status /= 0) error stop message
  write (*,"(DT(10))", IOSTAT=status, IOMSG=message) ndx
  if (status /= 0) error stop message

end subroutine ndx_test 

end program main