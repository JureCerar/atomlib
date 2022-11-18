# xyzio

The XYZ file format is a chemical file format. There is no formal standard and several variations exist. The file format is used in computational chemistry programs for importing and exporting geometries. The units are in Ångströms.
```
6
Two water molecules
O          1.26000       16.24000       16.79000
H          1.90000       16.61000       17.47000
H          1.77000       15.68000       16.13000
O         12.75000        0.53000        6.22000
H         13.37000        0.02000        6.80000
H         13.26000        1.20000        5.68000
```
Lines contain the following information (top to bottom):
- number of atoms
- comment string
- one line for each atom:
  - element name (optional)
  - x, y, and z coordinates (in Å)

For more information please see [reference](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html).

---------------------------------------------------------------------

## `xyz_open`
Open file for reading and checked. Frame offset table is constructed (byte position where frame starts). Move to selected frame with [fseek](https://gcc.gnu.org/onlinedocs/gfortran/FSEEK.html#FSEEK). If the function fails stat is a positive number, otherwise it is set to zero. If stat is not present end error occurs error stop is raised.
```fortran
call xyz_open (file, unit, nframes, offset, stat, errmsg)
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

## `xyz_close`
Close opened file.
```fortran
call xyz_close (unit, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `xyz_read_header`
Read file header in XYZ format.
```fortran
call xyz_read_header (unit, comment, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`character(*), intent(OUT) :: comment`
> Comment string.

`integer, intent(OUT) :: natoms`
> Number of atoms.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `xyz_read_data`
Read data in XYZ format.
```fortran
call xyz_read_data (unit, natoms, name, coor, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`character(*), intent(OUT) :: name`
> Atom name.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `xyz_read_coor`
Read ONLY coordinate data in XYZ format.
```fortran
call xyz_read_coor (unit, natoms, coor, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `xyz_skip_data`
Read through file (skip) in XYZ format but do not store any data.
```fortran
call xyz_skip_data (unit, natoms, stat, errmsg)
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

## `xyz_write`
Write data in XYZ format.
```fortran
call xyz_write (unit, comment, natoms, name, coor, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`character(*), intent(OUT) :: comment`
> Comment string.

`integer, intent(IN) :: natoms`
> Number of atoms.

`character(*), intent(IN) :: name`
> Atom name.

`real, dimension(DIM,natoms), intent(IN) :: coor`
> Atom cartesian coordinates in Å.

`real, dimension(DIM,natoms), intent(IN), optional :: vel`
> Atom velocity in nm/ps.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## Example
Example of reading XYZ file.
```fortran
program main
  use atomlib_xyzio
  use iso_fortran_env, only: INT64
  implicit none
  integer, parameter :: DIM = 3
  character(6), allocatable :: name(:)
  real, allocatable :: coor(:,:)
  integer(INT64), allocatable :: offset(:)
  character(128) :: comment
  integer :: unit, natoms, nframes, status

  call xyz_open ("conf.xyz", unit, nframes, offset, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xyz_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  allocate (name(natoms), coor(DIM,natoms), STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xyz_read_data (unit, natoms, name, coor, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xyz_close (unit, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message  

end program main
```