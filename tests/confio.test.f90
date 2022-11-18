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
  use atomlib_confio
  implicit none
  real, parameter :: DELTA = 0.001
  integer, parameter :: NUM = 10
  type(conf_t) :: conf
  character(128) :: arg, out, message
  integer :: i, status

  ! First check if all internal functions work properly.
  call conf_test ()

  ! Check for each file type
  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    write (*,*) "Processing file: '", trim(arg), "'"

    ! Load file and check values
    call conf%load(arg)
    if (conf%natoms /= 44) error stop

    select case (extension(arg))
    case ("xyz")
      if (conf%atomn(1) /= "N") error stop
      if (abs(conf%coor(1,1) - 0.630) > DELTA) error stop
    case ("gro")
      if (conf%atomi(1) /= 1) error stop
      if (conf%atomn(1) /= "N") error stop
      if (conf%resi(1) /= 1) error stop
      if (conf%resn(1) /= "ALA") error stop
      if (abs(conf%coor(1,1) - 0.630) > DELTA) error stop
      if (abs(conf%box(1,1) - 2.000) > DELTA) error stop
    case ("pdb")
      if (conf%atomi(1) /= 1) error stop
      if (conf%atomn(1) /= "N") error stop
      if (conf%resi(1) /= 1) error stop
      if (conf%resn(1) /= "ALA") error stop
      if (conf%chain(1) /= "A") error stop
      if (conf%element(1) /= "N") error stop
      if (abs(conf%bfact(1) - 0.000) > DELTA) error stop
      if (conf%charge(1) /= -1) error stop
      if (abs(conf%coor(1,1) - 0.630) > DELTA) error stop
      if (abs(conf%box(1,1) - 2.000) > DELTA) error stop
    case ("crd", "cor")
      if (conf%atomi(1) /= 1) error stop
      if (conf%atomn(1) /= "N") error stop
      if (conf%resi(1) /= 1) error stop
      if (conf%resn(1) /= "ALA") error stop
      if (conf%chain(1) /= "A") error stop
      if (abs(conf%bfact(1) - 1.000) > DELTA) error stop
      if (abs(conf%pcharge(1) - 0.000) > DELTA) error stop
      if (abs(conf%coor(1,1) - 0.630) > DELTA) error stop
    case default
      error stop "Unknown extension"
    end select

    ! Read file w/ error handling
    call conf%load(arg, STAT=status, ERRMSG=message)
    if (status /= 0) error stop
    
    ! Generate output file name
    out = "out." // extension(arg)
    
    ! Write file
    call conf%save (out)
    
    ! Write file w/ error handling
    call conf%save (out, STAT=status, ERRMSG=message)
    if (status /= 0) error stop

  end do

contains

! Test all internal functions of conf_t
subroutine conf_test ()
  implicit none
  integer, parameter :: NP = 10
  type(conf_t) :: conf, copy
  character(128) :: message
  integer :: i, status

  ! Allocate reallocate and check size
  call conf%allocate (NP)
  call conf%reallocate (2*NP)
  call conf%reallocate (NP/2)
  call conf%deallocate ()

  ! Allocate reallocate and check size
  call conf%allocate (NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call conf%reallocate (2*NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call conf%reallocate (NP/2, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call conf%deallocate (STAT=status, ERRMSG=message)
  if (status /= 0) error stop

  ! Construct a dummy configuration
  call conf%allocate(NP)
  conf%atomi   = [(i, i=1,NP)]
  conf%atomn   = "Xxx"
  conf%resi    = [(i, i=1,NP)]
  conf%resn    = "UNK"
  conf%element = "Xx"
  conf%chain   = "A"
  conf%bfact   = [(real(i), i=1,NP)]
  conf%pcharge = [(real(i), i=1,NP)]
  conf%charge  = [(i, i=1,NP)]
  conf%coor    = reshape([(real(i), i=1,3*NP)], [3,NP])
  conf%box     = reshape([(real(i), i=1,3*3)], [3,3])

  ! Construct a dummy configuration
  conf = conf_t(NP, ATOMI=[(i, i=1,NP)], ATOMN=[("Xxx", i=1,NP)], RESI=[(i, i=1,NP)], &
  & RESN=[("UNK", i=1,NP)], CHAIN=[("A", i=1,NP)], ELEMENT=[("X", i=1,NP)], BFACT=[(1.0, i=1,NP)], &
  & PCHARGE=[(1.0, i=1,NP)], CHARGE=[(1, i=1,NP)], COOR=reshape([(real(i),i=1,3*NP)],[3,NP]), &
  & BOX=reshape([(real(i),i=1,3*3)],[3,3]))

  ! Try copy and assign. All values should be copied and offset by NP.
  copy = conf + conf
  if (copy%atomi(1) /= copy%atomi(NP+1)) error stop
  if (copy%atomn(1) /= copy%atomn(NP+1)) error stop
  if (copy%resi(1) /= copy%resi(NP+1)) error stop
  if (copy%resn(1) /= copy%resn(NP+1)) error stop
  if (copy%element(1) /= copy%element(NP+1)) error stop
  if (copy%chain(1) /= copy%chain(NP+1)) error stop
  if (copy%bfact(1) /= copy%bfact(NP+1)) error stop
  if (copy%pcharge(1) /= copy%pcharge(NP+1)) error stop
  if (copy%charge(1) /= copy%charge(NP+1)) error stop
  if (copy%coor(1,1) /= copy%coor(1,NP+1)) error stop

  ! Try derived type write functions
  write (*, *, IOSTAT=status, IOMSG=message) conf
  if (status /= 0) error stop message
  write (*,"(DT)", IOSTAT=status, IOMSG=message) conf
  if (status /= 0) error stop message
  write (*,"(DT(17,8))", IOSTAT=status, IOMSG=message) conf
  if (status /= 0) error stop message

end subroutine conf_test 

! Util: Returns extension of file (always in lower case).
function extension (file) result (out)
  implicit none
  character(:), allocatable :: out
  character(*), intent(in) :: file
  integer, parameter :: offset = ichar("a") - ichar("A")
  integer :: i

  i = index(file, ".", BACK=.true.)
  if (i == 0) then
    out = ""
  else
    out = file(i+1:len_trim(file))
  end if

  do i = 1, len_trim(out)
    select case (out(i:i))
    case ("A" : "Z")
      out(i:i) = char(ichar(out(i:i)) - offset)
    end select
  end do

end function extension

end program main