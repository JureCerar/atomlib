# crdio

The CHARMM Cartesian Coordinate format (cor/crd) is standard format of CHARMM molecular modeling suit. This library only support non-binary file formats. The units are in Ångströms.

```
* TWO WATER MOLECULES
* BOX=  18.20600  18.20600  18.20600
* 
    6
    1    1 WAT  OW1    1.26000  16.24000  16.79000 SOLV A      1.00000
    2    1 WAT  HW2    1.90000  16.61000  17.47000 SOLV A      1.00000
    3    1 WAT  HW3    1.77000  15.68000  16.13000 SOLV A      1.00000
    4    2 WAT  OW1   12.75000   0.53000   6.22000 SOLV A      1.00000
    5    2 WAT  HW2   13.37000   0.02000   6.80000 SOLV A      1.00000
    6    2 WAT  HW3   13.26000   1.20000   5.68000 SOLV A      1.00000
```

Lines contain the following information (top to bottom):
- lines starting with `*` are title lines.
- number of atoms
- one line for each atom:
  - atom index,
  - residue index,
  - residue name,
  - atom name,
  - x, y, and z coordinates (in Å)
  - segment id,
  - residue id,
  - and weight.

Fortran string for less than 100000 atoms and PSF IDs with less than five characters `(i5,i5,x,a4,x,a4,3(f10.5),x,a4,x,a4,f10.5)`, otherwise `(i10,i10,2x,a8,2x,a8,3(f20.10),2x,a8,2x,a8,f20.10)`.

---------------------------------------------------------------------

## `crd_open`
Open file for reading and checked. Frame offset table is constructed (byte position where frame starts). Move to selected frame with [fseek](https://gcc.gnu.org/onlinedocs/gfortran/FSEEK.html#FSEEK). If the function fails stat is a positive number, otherwise it is set to zero. If stat is not present end error occurs error stop is raised.
```fortran
call crd_open (file, unit, nframes, offset, stat, errmsg)
```
#### *Parameters*
`character(*), intent(IN) :: file`
> Input file name.

`integer, intent(OUT) :: unit`
> Fortran file unit.

`integer, intent(OUT) :: nframes`
> Number of frames in file.

`integer(INT64), allocatable, dimension(:), intent(INOUT) :: offset`
> Frame offset table. Array size same as nframes.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `crd_close`
Close opened file.
```fortran
call crd_close (unit, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `crd_read_header`
Read file header in CRD/COR format.
```fortran
call crd_read_header (unit, extended, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`logical, intent(OUT) :: extended`
> Is file in normal or extended format.

`integer, intent(OUT) :: natoms`
> Number of atoms.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `crd_read_data`
Read data in PDB format.
```fortran
call crd_read_data (unit, extended, natoms, atomi, atomn, resi, resn, & 
& segid, resic, wfact, coor, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`logical, intent(IN) :: extended`
> Is file in normal or extended format.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, dimension(natoms), intent(OUT) :: atomi, resi`
> Atom and residue index.

`character(*), dimension(natoms), intent(OUT) :: atomn, resn`
> Atom and residue name.

`character(*), intent(OUT) :: segid, resic`
> Segment identifier and alternative residue identifier

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`real, dimension(natoms), intent(OUT) :: wfact`
> Weight factor.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `crd_read_coor`
Read ONLY coordinate data in CRD/COR format.
```fortran
call crd_read_coor (unit, extended, natoms, coor, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`logical, intent(IN) :: extended`
> Is file in normal or extended format.

`integer, intent(IN) :: natoms`
> Number of atoms.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `crd_skip_data`
Read through file (skip) in CRD/COR format but do not store any data.
```fortran
call crd_skip_data (unit, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `crd_write`
Read data in CRD/COR format. Function will automaticaly write in expanded format when needed.
```fortran
call pdb_write (unit, natoms, atomi, atomn, resi, resn, segid, resic, & 
& wfact, coor, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, dimension(natoms), intent(OUT) :: atomi, resi`
> Atom and residue index.

`character(*), dimension(natoms), intent(OUT) :: atomn, resn`
> Atom and residue name.

`character(*), intent(OUT) :: segid, resic`
> Segment identifier and alternative residue identifier

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`real, dimension(natoms), intent(OUT) :: wfact`
> Weight factor.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## Example
Example of reading PDB file.
```fortran
program main
  use atomlib_crdio
  use iso_fortran_env, only: INT64
  implicit none
  integer, parameter :: DIM = 3
  integer(INT64), allocatable :: offset(:)
  integer, allocatable :: atomi(:), resi(:)
  character(6), allocatable :: atomn(:), resn(:), segid(:), resic(:)
  real, allocatable :: wfact(:), coor(:,:)
  logical :: extension
  real :: box(DIM,DIM)
  integer :: i, unit, nframes, natoms, status

  call crd_open (arg, unit, nframes, offset, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call crd_read_header (unit, extended, natoms, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  allocate (atomi(natoms), atomn(natoms), resi(natoms), resn(natoms), segid(natoms), resic(natoms), &
  & wfact(natoms), coor(DIM,natoms))

  call crd_read_data (unit, extended, natoms, atomi, atomn, resi, resn, & 
  & segid, resic, wfact, coor, STAT=status, ERRMSG=message) 
  if (status /= 0) error stop message

  call crd_close (unit, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message 

end program main
```