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
  use atomlib_frame
  use atomlib_trajio
  implicit none
  integer, parameter :: SEEK_SET = 0
  real, parameter :: DELTA = 0.001
  type(traj_t) :: traj
  character(128) :: arg, message
  integer :: i, status

  ! First check if all internal functions work properly.
  call frame_test ()
  call traj_test ()

  ! Check for each file type
  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    write (*,*) "Processing file: '" // trim(arg) //"'"

    ! Open trajectory and move through frames.
    call traj%open(arg)
    call traj%read_next()
    call traj%read_next(2)
    call traj%fseek(1, SEEK_SET)
    if (traj%ftell() /= 1) error stop
    call traj%close()

    ! Open trajectory and move through frames /w error check.
    call traj%open(arg, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call traj%read_next(STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call traj%read_next(2, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call traj%fseek(1, SEEK_SET, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call traj%close(STAT=status, ERRMSG=message)
    if (status /= 0) error stop message

    ! Load file and check values
    call traj%load(arg)
    if (traj%nframes /= 3) error stop
    if (traj%frame(1)%natoms /= 44) error stop
    if (abs(traj%frame(1)%coor(1,1) - 0.630) > DELTA) error stop

    ! Load file with stride and offset
    call traj%load(arg, FIRST=1, LAST=3, STRIDE=2)
    if (traj%nframes /= 2) error stop

    ! Load file w/ error handling
    call traj%load(arg, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
  end do

  ! Write file to output
  call traj%save("out.xyz")
  call traj%save("out.dcd")
  call traj%save("out.trr")
  call traj%save("out.xtc")

  ! Write file to output /w error handling.
  if (status /= 0) error stop message
  call traj%save("out.xyz", STAT=status, ERRMSG=message)
  if (status /= 0) error stop message
  call traj%save("out.dcd", STAT=status, ERRMSG=message)
  if (status /= 0) error stop message
  call traj%save("out.trr", STAT=status, ERRMSG=message)
  if (status /= 0) error stop message
  call traj%save("out.xtc", STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

contains

! Test all internal functions of frame_t object.
subroutine frame_test ()
  implicit none
  integer, parameter :: NP = 10
  type(frame_t) :: frame, copy
  character(128) :: message
  integer :: i, status

  ! Allocate reallocate and check size
  call frame%allocate (NP)
  call frame%reallocate (2*NP)
  call frame%reallocate (NP/2)
  call frame%deallocate ()

  ! Allocate reallocate and check size
  call frame%allocate (NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call frame%reallocate (2*NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call frame%reallocate (NP/2, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call frame%deallocate (STAT=status, ERRMSG=message)
  if (status /= 0) error stop

  ! Construct a dummy configuration
  call frame%allocate(NP)
  frame%time = 1.0
  frame%coor = reshape([(real(i), i=1,3*NP)], [3,NP])
  frame%box = reshape([(real(i), i=1,3*3)], [3,3])

  ! Construct a dummy configuration
  frame = frame_t(NP, TIME=1.0, COOR=reshape([(real(i), i=1,3*NP)], [3,NP]), &
  & BOX=reshape([(real(i), i=1,3*3)], [3,3]))

  ! Try copy and assign. All values should be copied and offset by NP.
  copy = frame + frame
  if (copy%natoms /= 2*NP) error stop
  if (copy%time /= frame%time) error stop
  if (any(copy%box /= frame%box)) error stop
  if (copy%coor(1,1) /= copy%coor(1,NP+1)) error stop

  ! Try derived type write functions
  write (*, *, IOSTAT=status, IOMSG=message) frame
  if (status /= 0) error stop message
  write (*,"(DT)", IOSTAT=status, IOMSG=message) frame
  if (status /= 0) error stop message
  write (*,"(DT(17,8))", IOSTAT=status, IOMSG=message) frame
  if (status /= 0) error stop message

end subroutine frame_test 

! Test all internal functions of traj_t object.
subroutine traj_test ()
  implicit none
  integer, parameter :: NP = 4, NX = 10
  type(traj_t) :: traj
  type(frame_t) :: frame
  character(128) :: message
  integer :: i, status

  ! Allocate reallocate and check size
  call traj%allocate (NP)
  call traj%reallocate (2*NP)
  call traj%reallocate (NP/2)
  call traj%deallocate ()

  ! Allocate reallocate and check size
  call traj%allocate (NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call traj%reallocate (2*NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call traj%reallocate (NP/2, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call traj%deallocate (STAT=status, ERRMSG=message)
  if (status /= 0) error stop

  ! Construct a dummy configuration.
  call frame%allocate(NX) 
  frame%time = 1.0
  frame%coor = reshape([(real(i), i=1,3*NX)], [3,NX])
  frame%box = reshape([(real(i), i=1,3*3)], [3,3])

  ! Construct a dummy configuration
  call traj%allocate(NP)
  traj%frame = frame

  ! Construct a dummy configuration
  traj = traj_t(NP, FRAME=[(frame, i=1,NP)])

  ! Try derived type write functions
  write (*, *, IOSTAT=status, IOMSG=message) traj
  if (status /= 0) error stop message
  write (*,"(DT)", IOSTAT=status, IOMSG=message) traj
  if (status /= 0) error stop message
  write (*,"(DT(17,8))", IOSTAT=status, IOMSG=message) traj
  if (status /= 0) error stop message

end subroutine traj_test 

end program main