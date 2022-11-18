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

module atomlib_pdbio
  implicit none
  private
  public :: pdb_open, pdb_close
  public :: pdb_get_natoms, pdb_read_data, pdb_read_coor, pdb_skip_data, pdb_write

  ! TODO:
  ! - [ ] Implement "ORIGXn" and "SCALEn" directives.
  ! - [ ] Unify parameter order in functions.

  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension
  integer, parameter :: ERRLEN = 128 ! Lenght of error message
  integer, parameter :: PDB_MAX_LEN = 80 ! Maximum length of string
  integer, parameter :: PDB_MAX_NUM = 100000 ! Max. values of residue and atom index

  ! Structure of PDB file:
  !          1         2         3         4         5         6         7         8
  ! 12345678901234567890123456789012345678901234567890123456789012345678901234567890
  ! ATOM  12345  N   HIS A   1      49.668  24.248  10.436  1.00 25.00           N
  ! 
  ! COLUMNS   FMT   VARIABLE     DEFINITION
  !  1 -  6   a6                 Record type: "ATOM" or "HETATM" 
  !  7 - 11   i5    atomi        Atom serial number.
  ! 12        x                  NOTE: Additional character if atomname is long. 
  ! 13 - 16   a4    atomn        Atom name. 
  ! 17        a     altloc       Alternate location indicator. 
  ! 18 - 20   a3    resn         Residue name.
  ! 21        x                  NOTE: Additional character if resname is long.
  ! 22        a     chain        Chain identifier. 
  ! 23 - 26   i4    resi         Residue sequence number 
  ! 27        a     resic        Code for insertions of residues
  ! 31 - 38   f8.3  x            Orthogonal coordinates for X in Angstroms. 
  ! 39 - 46   f8.3  y            Orthogonal coordinates for Y in Angstroms. 
  ! 47 - 54   f8.3  z            Orthogonal coordinates for Z in Angstroms. 
  ! 55 - 60   f6.2  qfact        Occupancy. 
  ! 61 - 66   f6.2  bfact        Temperature factor.
  ! 73 - 76   x                  Segment identifier (obsolete). NOT IMPLEMENTED.
  ! 77 - 78   a2    element      Element symbol. 
  ! 79 - 80   i,a   charge       Charge on the atom (int + sign).
  !
  ! See: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
  
  ! order:
  ! unit, natoms, atomi, atomn, resi, resn, chain, resic, altloc, qfact, bfact, element, charge, coor, box

contains

! Open and check PDB file. Construct frame offset table. 
subroutine pdb_open (file, unit, nframes, offset, stat, errmsg)
  use iso_fortran_env, only: INT64, IOSTAT_END
  implicit none
  character(*), intent(in) :: file ! Name of file.
  integer, intent(out) :: unit ! Fortran file handle.
  integer, intent(out) :: nframes ! Num. of frames in file
  integer(INT64), allocatable, intent(out) :: offset(:) ! Stream position offsets (read mode only).
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer(INT64) :: pos
  integer(INT64), allocatable :: temp(:)
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
      inquire (UNIT=unit, POS=pos, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

      ! Read header and if EOF then there are (probably) no more frames in trajectory.
      call pdb_get_natoms (unit, natoms, STAT=status, ERRMSG=message)
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
      call pdb_skip_data (unit, natoms, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    end do while

    ! Rewind to start of the file for actual processing.
    rewind (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Reallocate offset table to actual size.
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

end subroutine pdb_open

! Close PDB file.
subroutine pdb_close (unit, stat, errmsg)
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

end subroutine pdb_close

! Util: Count num. of atoms in current PDB frame i.e. count occurance of "ATOM" and "HETATOM" keywords.
! NOTE: File must be opened as 'formatted stream'.
subroutine pdb_get_natoms (unit, natoms, stat, errmsg)
  use iso_fortran_env, only: INT64
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(out) :: natoms ! Number of atoms.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer(INT64) :: pos
  character(6) :: keyword
  character(9) :: action, access
  integer :: status
  logical :: opened

  natoms = 0

  catch: block

    ! Check if file unit is opened for reading and as stream.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    else if (index(access, "STREAM") == 0) then
      message = "File must be opened as stream"
      status = 1
    end if
    if (status /= 0) exit catch

    inquire (UNIT=unit, POS=pos, IOMSG=message, IOSTAT=status)
    if (status /= 0) exit catch

    ! Count occurance of keywords "ATOM" and "HETATM" until "END" keyword is encountered.
    do while (.true.)
      read (unit, "(a6)", IOSTAT=status, IOMSG=message) keyword
      if (status /= 0) exit catch
      select case (trim(keyword))
      case ("ATOM", "HETATM")
        natoms = natoms + 1
      case ("END", "ENDMDL")
        exit ! end of frame
      end select
    end do
  
    ! Return to original position using dummy stream read i.e. read 4 bytes before
    ! * the selected position. If position is at beggining simply rewind file.  
    if (pos < 4) then
      rewind (unit, IOSTAT=status, IOMSG=message)
    else
      read (unit, *, POS=pos-4, IOSTAT=status, IOMSG=message)
    end if
    if (status /= 0) exit catch

   ! Check if number of atoms makes sense
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

end subroutine pdb_get_natoms

! Read data in PDB format.
subroutine pdb_read_data (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
& altloc, qfact, bfact, element, charge, coor, box, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms 
  integer, intent(out) :: atomi(natoms) ! Atom index
  character(*), intent(out) :: atomn(natoms) ! Atom name
  integer, intent(out) :: resi(natoms) ! Residue index
  character(*), intent(out) :: resn(natoms) ! Residue name
  character(*), intent(out) :: chain(natoms) ! Chain identifier
  character(*), intent(out) :: resic(natoms) ! Code for insertions of residues.
  character(*), intent(out) :: altloc(natoms) ! Alternate location indicator. 
  real, intent(out) :: qfact(natoms) ! Occupancy factor.
  real, intent(out) :: bfact(natoms) ! Temperature (beta) factor.
  character(*), intent(out) :: element(natoms) ! Element symbol
  integer, intent(out) :: charge(natoms) ! Elemental (formal) charge of atoms
  real, intent(out) :: coor(DIM,natoms) ! Cartesian coordinates in [A] units.
  real, intent(out) :: box(DIM,DIM) ! Simulation box side in [A] units. 
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(PDB_MAX_LEN) :: buffer
  character(9) :: action, keyword
  logical :: opened
  character :: sign
  integer :: i, j, status

  catch: block

    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch

    box = 0.0
    
    ! Read line and process it according to keyword. Atom data is indicated by 
    ! * lines starting with "ATOM" or "HETATOM" keywords. "CRYST1" holds info on
    ! * crystal cell which can also be simulation box. "END" indicates end of frame
    ! * "i" indicates current atom index.
    i = 1
    do while (.true.)
      read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
      if (status /= 0) exit catch

      ! Process keyword
      keyword = buffer(1:6)
      select case (trim(keyword))
      case ("ATOM", "HETATM")

        if (i > natoms) then
          message = "Invalid number of atoms"
          stat = 1
          exit catch
        end if

        ! Index:     0  11  12  16  17  20  21  22  26  27  30     54     66   76  78      80
        100 format (6x, i5, 1x, a4, a1, a3, 1x, a1, i4, a1, 3x, 3f8.3, 2f6.2, 10x, a2, i1, a1)
        read (buffer, 100, IOSTAT=status, IOMSG=message) atomi(i), atomn(i), altloc(i), resn(i), chain(i), &
        & resi(i), resic(i), coor(:,i), qfact(i), bfact(i), element(i), charge(i), sign
        if (status /= 0) exit catch

        ! Append additional (extra) character to atom name and account for charge sign. 
        atomn(i) = trim(trim(atomn(i)(2:)) // atomn(i)(1:1))      
        if (sign == "-") charge(i) = -charge(i)
        element(i) = adjustl(element(i)) 

        i = i + 1

      case ("CRYST1")
        ! TODO: I am too lazy to implement any other box definition. Just use cubic box... Please.
        read (buffer, "(6x,3f9.3)", IOSTAT=status, IOMSG=message) (box(j,j), j = 1, DIM)
        if (status /= 0) exit catch

      case ("END", "ENDMDL")
        exit

      end select
    end do ! while

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdb_read_data

! Read data in PDB format and only store coordinates and box.
subroutine pdb_read_coor (unit, natoms, coor, box, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms 
  real, intent(out) :: coor(DIM,natoms) ! Cartesian coordinates in [A] units.
  real, intent(out) :: box(DIM,DIM) ! Simulation box side in [A] units. 
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(PDB_MAX_LEN) :: buffer
  character(9) :: action, keyword
  logical :: opened
  integer :: i, j, status

  catch: block

    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch

    box = 0.0

    ! Read line and process it according to keyword. Atom data is indicated by 
    ! * lines starting with "ATOM" or "HETATOM" keywords. "CRYST1" holds info on
    ! * crystal cell which can also be simulation box. "END" indicates end of frame
    ! * "i" indicates current atom index.
    i = 1
    do while (.true.)
      read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
      if (status /= 0) exit catch

      ! Process keyword
      keyword = buffer(1:6)
      select case (trim(keyword))
      case ("ATOM", "HETATM")

        if (i > natoms) then
          message = "Invalid number of atoms"
          stat = 1
          exit catch
        end if

        read (buffer, "(30x,3f8.3)", IOSTAT=status, IOMSG=message) coor(:,i)
        if (status /= 0) exit catch

        i = i + 1

      case ("CRYST1")
        ! TODO: I am too lazy to implement any other box definition. Just use cubic box... Please.
        read (buffer, "(6x,3f9.3)", IOSTAT=status, IOMSG=message) (box(j,j), j = 1, DIM)
        if (status /= 0) exit catch

      case ("END", "ENDMDL")
        exit

      end select
    end do ! while

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdb_read_coor

! Read through PDB file but do not store any data.
subroutine pdb_skip_data (unit, natoms, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms 
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(9) :: action, keyword
  logical :: opened
  integer :: i,  status

  catch: block

    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Skip all lines but count number of "ATOM" and "HETATM" to check if 
    ! * everything is OK. "i" indicates current atom index.
    i = 0
    do while (.true.)
      read (unit, "(a6)", IOSTAT=status, IOMSG=message) keyword
      if (status /= 0) exit catch
      
      select case (trim(keyword))
      case ("ATOM", "HETATM")
        ! Check if index is below limit
        if (i > natoms) then
          message = "Invalid number of atoms"
          status = 1
          exit catch
        end if
        i = i + 1

      case ("END", "ENDMDL")
        exit

      end select
    end do

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdb_skip_data

! Write data in PDB format.
subroutine pdb_write (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
& altloc, qfact, bfact, element, charge, coor, box, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms 
  integer, intent(in) :: atomi(natoms) ! Atom index
  character(*), intent(in) :: atomn(natoms) ! Atom name
  integer, intent(in) :: resi(natoms) ! Residue index
  character(*), intent(in) :: resn(natoms) ! Residue name
  character(*), intent(in) :: chain(natoms) ! Chain identifier
  character(*), intent(in) :: resic(natoms) ! Code for insertions of residues.
  character(*), intent(in) :: altloc(natoms) ! Alternate location indicator. 
  real, intent(in) :: qfact(natoms) ! Occupancy factor.
  real, intent(in) :: bfact(natoms) ! Temperature (beta) factor.
  character(*), intent(in) :: element(natoms) ! Element symbol
  integer, intent(in) :: charge(natoms) ! Elemental (formal) charge of atoms
  real, intent(in) :: coor(DIM,natoms) ! Cartesian coordinates in [A] units.
  real, intent(in) :: box(DIM,DIM) ! Simulation box side in [A] units. 
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(6), parameter :: RECORD = "ATOM"
  character(8) :: resn_, atomn_, element_, charge_
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

    ! TODO: I am too lazy to implement any other box definition. Just use cubic box... Please.
    100 format ("CRYST1", 3f9.3, 3f7.2, x, a10, i3)
    write (unit, 100, IOSTAT=status, IOMSG=message) (box(i,i), i = 1, DIM), 90., 90., 90., "P 1", 1
    if (status /= 0) exit catch

    ! Write data frame. But first do some preprocessing. If atom name is longer than 4 characters
    ! * add last character to the beginning (additional character). Also generate elemental charge
    ! * in correct format i.e. "1+" and "1-"
    do i = 1, natoms
      ! Do some pre-processing
      atomn_ = atomn(i)(4:4) // atomn(i)(1:3)
      resn_ = adjustr(resn(i)(1:3))
      element_ = adjustr(element(i)(1:2))
      if (charge(i) /= 0) then
        write (charge_, "(i1,a1)") abs(charge(i)), merge("+", "-", charge(i) > 0)
      else
        charge_ = ""
      end if

      ! Index:     0  11  12  16  17  20  21  22  26  27  30     54     66   76  78  80
      200 format (a6, i5, 1x, a4, a1, a3, 1x, a1, i4, a1, 3x, 3f8.3, 2f6.2, 10x, a2, a2)
      write (unit, 200, IOSTAT=status, IOMSG=message) RECORD, mod(atomi(i), PDB_MAX_NUM), atomn_, altloc(i), &
      & resn(i), chain(i), mod(resi(i), PDB_MAX_NUM), resic(i), coor(:,i), qfact(i), bfact(i), element_, charge_
      if (status /= 0) exit catch

    end do ! for i

    write (unit, "(a,/,a)", IOSTAT=status, IOMSG=message) "TER", "END"
    if (status /= 0) exit catch

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine pdb_write

end module atomlib_pdbio
