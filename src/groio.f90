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

module atomlib_groio
  implicit none
  private
  public :: gro_open, gro_close, gro_read_header, gro_read_data, gro_read_coor, gro_skip_data, gro_write

  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension
  integer, parameter :: ERRLEN = 256 ! Length of error string
  integer, parameter :: GRO_MAX_NUM = 100000 ! Max. values of residue and atom index

  ! Structure of GRO file:
  ! See: https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#gro 
  !
  ! COMMENT
  ! NATOMS
  ! RESI RESN ATOMN ATOMI  X Y Z  Vx Vy Vz (velocities are optional)
  !   i5   a5    a5    i5  3f8.3     3f8.4
  ! ...
  ! SIDEx SIDEy SIDEz
  !  f9.5  f9.5  f9.5
  

contains

! Open and check GRO file. Construct frame offset table. 
subroutine gro_open (file, unit, nframes, offset, stat, errmsg)
  use iso_fortran_env, only: INT64, IOSTAT_END
  implicit none
  character(*), intent(in) :: file ! File name.
  integer, intent(out) :: unit ! Fortran file handle.
  integer, intent(out) :: nframes ! Num. of frames in file
  integer(INT64), allocatable, intent(out) :: offset(:) ! Stream position offsets.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer(INT64) :: pos
  integer(INT64), allocatable :: temp(:)
  character(64) :: comment
  integer :: status, natoms

  catch: block

    ! Open file as stream.
    open (FILE=file, NEWUNIT=unit,  ACCESS="stream", FORM="formatted", &
    &     STATUS="old", ACTION="read", IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Start at offset table of size 128. If needed, we will expand it as we go along.
    if (allocated(offset)) deallocate (offset, STAT=status, ERRMSG=message)
    allocate (offset(128), SOURCE=0_INT64, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    nframes = 0

    while: do while (.true.)
      ! Store current position
      inquire (UNIT=unit, POS=pos, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

      ! Read header and if EOF then there are (probably) no more frames in trajectory.
      call gro_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
      if (status == IOSTAT_END) exit while
      if (status /= 0) exit catch

      ! Expand offset table if needed
      if (nframes > size(offset)) then
        allocate (temp(2 * size(offset)), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        temp(1:nframes) = offset(1:nframes)
        call move_alloc (temp, offset)
      end if
      
      nframes = nframes + 1
      offset(nframes) = pos

      ! Skip remaining frame and only check for errors.
      call gro_skip_data (unit, natoms, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    end do while

    ! Rewind to start of the file for actual processing.
    rewind (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Allocate offset table to actual size.
    allocate (temp(nframes), STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    temp = offset(1:nframes)
    call move_alloc (temp, offset)

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine gro_open

! Close GRO file.
subroutine gro_close (unit, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  close (UNIT=unit, IOSTAT=status, IOMSG=message)
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message) 

end subroutine gro_close

! Read header in GRO format.
subroutine gro_read_header (unit, comment, natoms, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  character(*), intent(out) :: comment ! Comment line
  integer, intent(out) :: natoms ! Number of atoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(9) :: action
  logical :: opened
  integer :: status

  catch: block
    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch

    read (unit, "(a)", IOSTAT=status, IOMSG=message) comment
    if (status /= 0) exit catch    

    read (unit, *, IOSTAT=status, IOMSG=message) natoms
    if (status /= 0) exit catch 

    if (natoms < 0) then
      message = "Invalid number of atoms"
      status = 1
      exit catch
    end if

  end block catch
 
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)  

end subroutine gro_read_header 

! Read data in GRO format.
subroutine gro_read_data (unit, natoms, atomi, atomn, resi, resn, coor, box, vel, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms
  integer, intent(out) :: resi(natoms) ! Residue index
  character(*), intent(out) :: resn(natoms) ! Residue name
  integer, intent(out) :: atomi(natoms) ! Particle index
  character(*), intent(out) :: atomn(natoms) ! Particle name
  real, intent(out) :: coor(DIM,natoms) ! Cartesian coordinates in [nm] units.
  real, intent(out) :: box(DIM,DIM) ! Simulation box side in [nm] units.
  real, intent(out), optional :: vel(DIM,natoms) !  Velocities in [nm/ps] units. 
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(128) :: buffer
  character(9) :: action
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch
    
    ! If velocities are provided read everything, otherwise skip them.
    if (present(vel)) then
      do i = 1, natoms
        100 format (i5,2a5,i5,3f8.3,3f8.4)
        read (unit, 100, IOSTAT=status, IOMSG=message) resi(i), resn(i), atomn(i), atomi(i), coor(:,i), vel(:,i)
        if (status /= 0) exit catch
      end do

    else
      do i = 1, natoms
        200 format (i5,2a5,i5,3f8.3)
        read (unit, 200, IOSTAT=status, IOMSG=message) resi(i), resn(i), atomn(i), atomi(i), coor(:,i)
        if (status /= 0) exit catch
      end do

    end if

    ! Atom names are right-adjusted
    atomn(:) = adjustl(atomn)

    ! Read box dimensions and store them in buffer. First try to read as DIM vector (triclinic box).
    ! * If it fails, try to read then as DIM*DIM tensor (cubic box). If that fails then we have a problem. 
    read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
    if (status /= 0) exit catch

    read (buffer, *, IOSTAT=status) box(:,:) ! Triclinic box
    if (status /= 0) then
      box = 0.
      read (buffer, *, IOSTAT=status, IOMSG=message) (box(i,i), i = 1, DIM) ! Cubic box
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine gro_read_data

! Read data in GRO format and only store coordinates and box.
subroutine gro_read_coor (unit, natoms, coor, box, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms
  real, intent(out) :: coor(DIM,natoms) ! Cartesian coordinates in [nm] units.
  real, intent(out) :: box(DIM,DIM) ! Simulation box side in [nm] units.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(128) :: buffer
  character(9) :: action
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch
    
    do i = 1, natoms
      read (unit, "(20x,3f8.3)", IOSTAT=status, IOMSG=message) coor(:,i)
      if (status /= 0) exit catch
    end do

    ! Read box dimensions and store them in buffer. First try to read as DIM vector (triclinic box).
    ! * If it fails, try to read then as DIM*DIM tensor (cubic box). If that fails then we have a problem. 
    read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
    if (status /= 0) exit catch

    read (buffer, *, IOSTAT=status) box(:,:) ! Triclinic box
    if (status /= 0) then
      box = 0.
      read (buffer, *, IOSTAT=status, IOMSG=message) (box(i,i), i = 1, DIM) ! Cubic box
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine gro_read_coor

! Read through GRO file but do not store any data.
subroutine gro_skip_data (unit, natoms, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! File unit
  integer, intent(in) :: natoms ! Number of atoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(9) :: action
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch
    
    if (natoms < 0) then
      message = "Invalid number of atoms"
      status = 1
      exit catch
    end if

    ! Skip atom data and simulation box side (+1).
    do i = 1, natoms + 1
      read (unit, *, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch
    end do

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine gro_skip_data

! Write data in GRO format.
subroutine gro_write (unit, comment, natoms, atomi, atomn, resi, resn, coor, box, vel, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  character(*), intent(in) :: comment ! Comment line
  integer, intent(in) :: natoms ! Number of atoms
  integer, intent(in) :: atomi(natoms) ! Particle index
  character(*), intent(in) :: atomn(natoms) ! Particle name
  integer, intent(in) :: resi(natoms) ! Residue index
  character(*), intent(in) :: resn(natoms) ! Residue name
  real, intent(in) :: coor(DIM,natoms) ! Cartesian coordinates in [nm] units.
  real, intent(in) :: box(DIM,DIM) ! Simulation box side in [nm] units.
  real, intent(in), optional :: vel(DIM,natoms) !  Velocities in [nm/ps] units. 
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(9) :: action 
  logical :: opened
  integer :: i, status

  catch: block
  
    ! Check if file unit is opened for writing.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "WRITE") == 0) then
      message = "File unit not opened for writing"
      status = 1
    end if
    if (status /= 0) exit catch

    write (unit, "(a)", IOSTAT=status, IOMSG=message) trim(comment)
    if (status /= 0) exit catch    

    write (unit, "(2x,i0)", IOSTAT=status, IOMSG=message) natoms
    if (status /= 0) exit catch

    ! If velocities are provided write everything otherwise skip them. Maximum value of atom or residue
    ! * index is 100000, so wrap around if bigger. Additionaly names are right adjusted.
    if (present(vel)) then
      do i = 1, natoms
        write (unit, "(i5,2a5,i5,3f8.3,3f8.4)", IOSTAT=status, IOMSG=message) mod(resi(i), GRO_MAX_NUM), &
        & adjustl(resn(i)), adjustr(trim(atomn(i))), mod(atomi(i), GRO_MAX_NUM), coor(:,i), vel(:,i)
        if (status /= 0) exit catch
      end do

    else
      do i = 1, natoms
        write (unit, "(i5,2a5,i5,3f8.3)", IOSTAT=status, IOMSG=message) mod(resi(i), GRO_MAX_NUM), &
        & adjustl(resn(i)), adjustr(trim(atomn(i))), mod(atomi(i), GRO_MAX_NUM), coor(:,i)
        if (status /= 0) exit catch
      end do

    end if 

    ! If only diagonal values of simulation box are present write as DIM vector (triclinic box), otherwise
    ! * write as DIM*DIM tensor (cubic box).
    if (all(box /= 0.0)) then
      write (unit, "(*(x,f9.5))", IOSTAT=status, IOMSG=message) box(:,:)
      if (status /= 0) exit catch
    else
      write (unit, "(*(x,f9.5))", IOSTAT=status, IOMSG=message) (box(i,i), i = 1, DIM)
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine gro_write

end module atomlib_groio
