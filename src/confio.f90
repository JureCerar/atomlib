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

module atomlib_confio
  use atomlib_xyzio
  use atomlib_groio
  use atomlib_pdbio
  use atomlib_crdio
  implicit none
  private
  public :: conf_t

  integer, parameter :: DIM = 3
  integer, parameter :: MAXLEN = 8 
  integer, parameter :: ERRLEN = 128

  type :: conf_t
    integer, allocatable :: atomi(:) ! Atom index
    character(MAXLEN), allocatable :: atomn(:) ! Atom name
    integer, allocatable :: resi(:) ! Residue index
    character(MAXLEN), allocatable :: resn(:) ! Residue name
    character(MAXLEN), allocatable :: chain(:) ! Chain name
    character(MAXLEN), allocatable :: element(:) ! Element name
    real, allocatable :: bfact(:) ! Beta-factor
    real, allocatable :: pcharge(:) ! Partial-charge
    integer, allocatable :: charge(:) ! Formal-charge
    real, allocatable :: coor(:,:) ! Atom cartesian coordinates
    real :: box(DIM,DIM) = 0.0 ! Simulation box vector
    integer :: natoms = 0 ! Num. of atoms 
  contains
    ! Allocation functions
    procedure :: allocate => conf_allocate
    procedure :: reallocate => conf_reallocate
    procedure :: deallocate => conf_deallocate
    procedure :: allocated => conf_allocated
    ! Assigment/copy functions
    procedure, private :: conf_append, conf_assign
    generic :: operator(+) => conf_append
    generic :: assignment(=) => conf_assign
    ! Formatted I/O
    procedure, private :: conf_write_formatted
    generic :: write(formatted) => conf_write_formatted
    ! File I/O
    procedure :: load => conf_load
    procedure :: save => conf_save
  end type conf_t

  ! Constructor interface
  interface conf_t
    procedure :: conf_constructor
  end interface conf_t

contains

! Util: Returns extension of file (always in upper case).
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
    case ("a" : "z")
      out(i:i) = char(ichar(out(i:i)) - offset)
    end select
  end do

end function extension

! Constructor for conf_t dferived type.
type(conf_t) function conf_constructor (natoms, atomi, atomn, resi, resn, chain, element, &
& bfact, pcharge, charge, coor, box) result (this)
  implicit none
  integer, intent(in) :: natoms
  integer, intent(in), optional :: atomi(natoms)
  character(*), intent(in), optional :: atomn(natoms)
  integer, intent(in), optional :: resi(natoms)
  character(*), intent(in), optional :: resn(natoms)
  character(*), intent(in), optional :: chain(natoms)
  character(*), intent(in), optional :: element(natoms)
  real, intent(in), optional :: bfact(natoms)
  real, intent(in), optional :: pcharge(natoms)
  integer, intent(in), optional :: charge(natoms)
  real, intent(in), optional :: coor(DIM,natoms)
  real, intent(in), optional :: box(DIM,DIM)

  call this%allocate(natoms)
  if (present(atomi)) this%atomi = atomi
  if (present(atomn)) this%atomn = atomn
  if (present(resi)) this%resi = resi
  if (present(resn)) this%resn = resn
  if (present(element)) this%element = element
  if (present(chain)) this%chain = chain
  if (present(bfact)) this%bfact = bfact
  if (present(pcharge)) this%pcharge = pcharge
  if (present(charge)) this%charge = charge
  if (present(coor)) this%coor = coor
  if (present(box)) this%box = box

end function conf_constructor

! Allocate and initialize object data for specified number of atoms.
subroutine conf_allocate (this, natoms, stat, errmsg)
  implicit none
  class(conf_t), intent(inout) :: this
  integer, intent(in) :: natoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(MAXLEN), parameter :: empty = ""
  character(ERRLEN) :: message
  integer :: status

  catch: block

    if (natoms <= 0) then
      message = "Invalid number of particles"
      status = 1
      exit catch
    end if

    this%natoms = natoms

    if (allocated(this%atomi)) deallocate(this%atomi, STAT=status, ERRMSG=message)
    allocate (this%atomi(natoms), SOURCE=0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%atomn)) deallocate(this%atomn, STAT=status, ERRMSG=message)
    allocate (this%atomn(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%resi)) deallocate(this%resi, STAT=status, ERRMSG=message)
    allocate (this%resi(natoms), SOURCE=0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%resn)) deallocate(this%resn, STAT=status, ERRMSG=message)
    allocate (this%resn(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%chain)) deallocate(this%chain, STAT=status, ERRMSG=message)
    allocate (this%chain(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    
    if (allocated(this%element)) deallocate(this%element, STAT=status, ERRMSG=message)
    allocate (this%element(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%bfact)) deallocate(this%bfact, STAT=status, ERRMSG=message)
    allocate (this%bfact(natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%pcharge)) deallocate(this%pcharge, STAT=status, ERRMSG=message)
    allocate (this%pcharge(natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%charge)) deallocate(this%charge, STAT=status, ERRMSG=message)
    allocate (this%charge(natoms), SOURCE=0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%coor)) deallocate(this%coor, STAT=status, ERRMSG=message)
    allocate (this%coor(DIM, natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine conf_allocate

! Reallocate object data for specified number of atoms (bigger or smaller) in configuration.
subroutine conf_reallocate (this, natoms, stat, errmsg)
  implicit none
  class(conf_t), intent(inout) :: this
  integer, intent(in) :: natoms
  integer, intent(out), optional :: stat
  character(*), intent(out), optional :: errmsg
  character(MAXLEN), parameter :: empty = ""
  real, allocatable :: ftmp(:), f2tmp(:,:)
  integer, allocatable :: itmp(:)
  character(MAXLEN), allocatable :: ctmp(:)
  character(ERRLEN) :: message
  integer :: i, status
  
  catch: block 

    if (natoms <= 0) then
      message = "Invalid number of particles"
      status = 1
      exit catch
    end if

    this%natoms = natoms

    if (allocated(this%atomi)) then
      allocate (itmp(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%atomi))
      itmp(:i) = this%atomi(:i)
      call move_alloc (itmp, this%atomi)
    else
      allocate (this%atomi(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%atomn)) then
      allocate (ctmp(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%atomn))
      ctmp(:i) = this%atomn(:i)
      call move_alloc (ctmp, this%atomn)
    else
      allocate (this%atomn(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%resi)) then
      allocate (itmp(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%resi))
      itmp(:i) = this%resi(:i)
      call move_alloc (itmp, this%resi)
    else
      allocate (this%resi(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%resn)) then
      allocate (ctmp(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%resn))
      ctmp(:i) = this%resn(:i)
      call move_alloc (ctmp, this%resn)
    else
      allocate (this%resn(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if
    
    if (allocated(this%chain)) then
      allocate (ctmp(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%chain))
      ctmp(:i) = this%chain(:i)
      call move_alloc (ctmp, this%chain)
    else
      allocate (this%chain(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%element)) then
      allocate (ctmp(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%element))
      ctmp(:i) = this%element(:i)
      call move_alloc (ctmp, this%element)
    else
      allocate (this%element(natoms), SOURCE=empty, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%bfact)) then
      allocate (ftmp(natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%bfact))
      ftmp(:i) = this%bfact(:i)
      call move_alloc (ftmp, this%bfact)
    else
      allocate (this%bfact(natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%pcharge)) then
      allocate (ftmp(natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%pcharge))
      ftmp(:i) = this%pcharge(:i)
      call move_alloc (ftmp, this%pcharge)
    else
      allocate (this%pcharge(natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%charge)) then
      allocate (itmp(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%charge))
      itmp(:i) = this%charge(:i)
      call move_alloc (itmp, this%charge)
    else
      allocate (this%charge(natoms), SOURCE=0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

    if (allocated(this%coor)) then
      allocate (f2tmp(DIM, natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(natoms, size(this%coor, DIM=2))
      f2tmp(:,:i) = this%coor(:,:i)
      call move_alloc (f2tmp, this%coor)
    else
      allocate (this%coor(DIM, natoms), SOURCE=0.0, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if
    
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine conf_reallocate

! Deallocate object data.
subroutine conf_deallocate (this, stat, errmsg)
  implicit none
  class(conf_t), intent(inout) :: this
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  catch: block

    status = 0

    this%natoms = 0

    if (allocated(this%atomi)) deallocate(this%atomi, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%atomn)) deallocate(this%atomn, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%resi)) deallocate(this%resi, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%resn)) deallocate(this%resn, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    
    if (allocated(this%chain)) deallocate(this%chain, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%element)) deallocate(this%element, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%bfact)) deallocate(this%bfact, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%pcharge)) deallocate(this%pcharge, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%charge)) deallocate(this%charge, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    if (allocated(this%coor)) deallocate(this%coor, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine conf_deallocate

! Check if object is (properly) allocated.
function conf_allocated (this) result (out)
  implicit none  
  logical :: out
  class(conf_t), intent(in) :: this
  
  out = .true.
  out = out .and. allocated(this%atomi)
  out = out .and. allocated(this%atomn)
  out = out .and. allocated(this%resi)
  out = out .and. allocated(this%resn)
  out = out .and. allocated(this%chain)
  out = out .and. allocated(this%element)
  out = out .and. allocated(this%bfact)
  out = out .and. allocated(this%pcharge)
  out = out .and. allocated(this%charge)
  out = out .and. allocated(this%coor)

end function conf_allocated

! Append configuration to existing configuration. Both objects must be allocated.
type(conf_t) function conf_append (this, other) result (result)
  implicit none
  class(conf_t), intent(in) :: this
  type(conf_t), intent(in) :: other
  integer :: offset, natoms

  if (.not. (this%allocated() .and. other%allocated())) then
    error stop "Objects not allocated"
  end if 
  
  offset = size(this%atomi)
  natoms = offset + size(other%atomi)

  call result%allocate(natoms)

  result%atomi(:offset) = this%atomi
  result%atomn(:offset) = this%atomn
  result%resi(:offset) = this%resi
  result%resn(:offset) = this%resn
  result%chain(:offset) = this%chain
  result%element(:offset) = this%element
  result%bfact(:offset) = this%bfact
  result%pcharge(:offset) = this%pcharge
  result%charge(:offset) = this%charge
  result%coor(:,:offset) = this%coor
  result%box = this%box

  result%atomi(offset+1:) = other%atomi
  result%atomn(offset+1:) = other%atomn
  result%resi(offset+1:) = other%resi
  result%resn(offset+1:) = other%resn
  result%chain(offset+1:) = other%chain
  result%element(offset+1:) = other%element
  result%bfact(offset+1:) = other%bfact
  result%pcharge(offset+1:) = other%pcharge
  result%charge(offset+1:) = other%charge
  result%coor(:,offset+1:) = other%coor

  ! result%box = other%box

end function conf_append

! Assign value to configuration. Dump data if already allocated.
subroutine conf_assign (this, other)
  implicit none
  class(conf_t), intent(inout) :: this
  type(conf_t), intent(in) :: other
  integer :: natoms
  
  if (.not. other%allocated()) then
    error stop "Objects not allocated"
  end if 

  natoms = size(other%atomi)
  call this%allocate (natoms)

  this%atomi = other%atomi
  this%atomn = other%atomn
  this%resi = other%resi
  this%resn = other%resn
  this%chain = other%chain
  this%element = other%element
  this%bfact = other%bfact
  this%pcharge = other%pcharge
  this%charge = other%charge
  this%coor = other%coor
  this%box = other%box

end subroutine conf_assign

! Formatted write for configurations.
subroutine conf_write_formatted (this, unit, iotype, v_list, iostat, iomsg)
  implicit none
  class(conf_t), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype 
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  character(128) :: fmt, dt
  integer :: i

  catch: block

    if (.not. this%allocated()) then
      iomsg = "Object not allocated"
      iostat = 1
      exit catch
    end if 

    if (iotype == "LISTDIRECTED") then
      ! Number of atoms
      write (unit, "(i0,/)", IOSTAT=iostat, IOMSG=iomsg) this%natoms
      if (iostat /= 0) exit catch

      ! Data format:
      !    i8,  x,a8,  x,i8, x,a8,  x,a2,    x,a2, x,f8.4,  x,f8.4, SP,i3,S, 3(f8.4)
      ! atomi, atomn,  resi, resn, chain, element,  bfact, pcharge,  charge,    coor
      fmt = "(i8,x,a8,x,i8,x,a8,x,a2,x,a2,x,f8.4,x,f8.4,SP,i3,S,3(f8.4)/)"
      do i = 1, size(this%atomi)
        write (unit, fmt, IOSTAT=iostat, IOMSG=iomsg) this%atomi(i), this%atomn(i), this%resi(i), this%resn(i), this%chain(i), & 
        & this%element(i), this%bfact(i), this%pcharge(i), this%charge(i), this%coor(:,i)
        if (iostat /= 0) exit catch
      end do

      ! Simulation box
      write (unit, "(*(x,f8.4),/,/)", IOSTAT=iostat, IOMSG=iomsg) this%box
      if (iostat /= 0) exit catch

    else if (iotype == "DT") then
      if (size(v_list) == 0) then
        dt = "f8.4"
      else if (size(v_list) == 2) then
        write (dt,"(a,i0,a,i0)") "f", v_list(1), ".", v_list(2)
      else
        iostat = 1
        iomsg = "Integer-list for DT descriptor must contain two integers"
        exit catch
      end if

      ! Number of atoms
      write (unit, "(i0,/)", IOSTAT=iostat, IOMSG=iomsg) this%natoms
      if (iostat /= 0) exit catch

      ! Data box in user defined format.
      fmt = "(i8,x,a8,x,i8,x,a8,x,a2,x,a2,x,"//trim(dt)//","//trim(dt)//",SP,i3,S,3("//trim(dt)//"),/)"
      do i = 1, size(this%atomi)
        write (unit, fmt, IOSTAT=iostat, IOMSG=iomsg) this%atomi(i), this%atomn(i), this%resi(i), this%resn(i), this%chain(i), & 
        & this%element(i), this%bfact(i), this%pcharge(i), this%charge(i), this%coor(:,i)
        if (iostat /= 0) exit catch
      end do

      ! Simulation box in user defined format.
      fmt = "(*(x,"//trim(dt)//"),/,/)"
      write (unit, fmt, IOSTAT=iostat, IOMSG=iomsg) this%box
      if (iostat /= 0) exit catch
      
    else
      iomsg = "Unsupported iotype"
      iostat = 1
      exit catch

    end if

  end block catch

end subroutine conf_write_formatted

! Save (write) configuration to file. Format is based on file extension.
subroutine conf_save (this, file, stat, errmsg)
  implicit none
  class(conf_t), intent(in) :: this
  character(*), intent(in) :: file
  integer, intent(out), optional :: stat
  character(*), intent(out), optional :: errmsg
  character(ERRLEN) :: message
  integer :: i, unit, natoms, status

  catch: block

    if (.not. this%allocated()) then
      message = "Object not allocated"
      status = 1
      exit catch
    end if

    natoms = size(this%atomi)

    select case (extension(file))
    case ("XYZ")
      ! NOTE: XYZ default units are Angstroms.
      xyz: block
        character(80) :: comment = ""
        open (FILE=file, NEWUNIT=unit, STATUS="unknown", ACTION="write", IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        call xyz_write (unit, comment, natoms, this%atomn, this%coor * 10, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
      end block xyz

    case ("GRO")
      gro: block
        character(80) :: comment = ""
        open (FILE=file, NEWUNIT=unit, STATUS="unknown", ACTION="write", IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        call gro_write (unit, comment, natoms, this%atomi, this%atomn, this%resi, & 
        & this%resn, this%coor, this%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch   
      end block gro

    case ("PDB")
      ! NOTE: PDB default units are Angstroms.
      pdb: block
        open (FILE=file, NEWUNIT=unit, STATUS="unknown", ACTION="write", IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        call pdb_write (unit, natoms, this%atomi, this%atomn, this%resi, this%resn, this%chain,  &
        & [("", i = 1, natoms)], [("", i = 1, natoms)], [(0.0, i = 1, natoms)], this%bfact, this%element, &
        & this%charge, this%coor * 10, this%box * 10, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch         
      end block pdb

    case ("CRD", "COR")
      ! NOTE: CRD/COR default units are Angstroms.
      crd: block
        open (FILE=file, NEWUNIT=unit, STATUS="unknown", ACTION="write", IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        call crd_write (unit, natoms, this%atomi, this%atomn, this%resi, this%resn, this%chain, &
        & [("", i = 1, natoms)], this%bfact, this%coor * 10, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch     
      end block crd

    case default 
      message = "Unsuppoted file type: '" // extension(file) // "'"
      status = 1
      exit catch

    end select

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine conf_save

! Load (read) configuration file.
subroutine conf_load (this, file, stat, errmsg)
  implicit none
  class(conf_t), intent(inout) :: this
  character(*), intent(in) :: file
  integer, intent(out), optional :: stat
  character(*), intent(out), optional :: errmsg
  character(ERRLEN) :: message
  integer :: unit, natoms, status

  catch: block
  
    ! Read input according to data format.
    select case (extension(file))
    case ("XYZ")
      ! NOTE: .XYZ default units are Angstroms.
      xyz: block
        character(80) :: comment
        open (FILE=file, NEWUNIT=unit,  ACCESS="stream", FORM="formatted", &
        & STATUS="old", ACTION="read", IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        call xyz_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call this%allocate(natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call xyz_read_data (unit, natoms, this%atomn, this%coor, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        this%coor = this%coor * 0.1 ! to [nm].
        this%box = 0.0 ! No box information
      end block xyz

    case ("GRO")
      gro: block
        character(80) :: comment
        open (FILE=file, NEWUNIT=unit,  ACCESS="stream", FORM="formatted", &
        & STATUS="old", ACTION="read", IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        call gro_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call this%allocate(natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call gro_read_data (unit, natoms, this%atomi, this%atomn, this%resi, this%resn, &
        & this%coor, this%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
      end block gro

    case ("PDB")
      ! NOTE: .PDB default units are Angstroms.
      pdb: block
        character(8), allocatable :: dmy1(:), dmy2(:)
        real, allocatable :: dmy3(:)
        open (FILE=file, NEWUNIT=unit,  ACCESS="stream", FORM="formatted", &
        & STATUS="old", ACTION="read", IOSTAT=status, IOMSG=message)
        call pdb_get_natoms (unit, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        allocate (dmy1(natoms), dmy2(natoms), dmy3(natoms), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call this%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch 
        call pdb_read_data (unit, natoms, this%atomi, this%atomn, this%resi, this%resn, this%chain, dmy1, dmy2, dmy3, &
        & this%bfact, this%element, this%charge, this%coor, this%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        this%coor = this%coor * 0.1 ! to [nm].
        this%box = this%box * 0.1 ! to [nm].
      end block pdb

    case ("CRD", "COR")
      ! NOTE: .CRD/.COR default units are Angstroms.
      crd: block
        character(8), allocatable :: dmy(:)
        logical :: extended
        open (FILE=file, NEWUNIT=unit,  ACCESS="stream", FORM="formatted", &
        & STATUS="old", ACTION="read", IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        call crd_read_header (unit, extended, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call this%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        allocate (dmy(natoms), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call crd_read_data (unit, extended, natoms, this%atomi, this%atomn, this%resi, this%resn, & 
        & dmy, this%chain, this%bfact, this%coor, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        close (unit, IOSTAT=status, IOMSG=message)
        if (status /= 0) exit catch
        this%coor = this%coor * 0.1 ! to [nm].
        this%box = 0.0 ! No box information
      end block crd

    case default 
      message = "Unsuppoted file type: '" // extension(file) // "'"
      status = 1
      exit catch

    end select

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine conf_load

end module atomlib_confio
