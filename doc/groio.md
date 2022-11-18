# groio
Files with the gro file extension contain a molecular structure in Gromos87 format. gro files can be used as trajectory by simply concatenating files. 
```
Two water molecules, t= 0.0
    6
    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180
    2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734
    2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257
    2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244
   1.82060   1.82060   1.82060
```
Lines contain the following information (top to bottom):
- title string (free format string, optional time in ps after 't=') 
- number of atoms (free format)
- one line for each atom (fixed format):
  -  residue number
  -  residue name
  -  atom name
  -  atom number
  -  position (in nm, x y z in 3 columns)
  -  velocity (in nm/ps (or km/s), x y z in 3 columns) 
- box vectors (free format, space separated reals), values: vx, vy, bz

Fortran format: `(i5,2a5,i5,3f8.3,3f8.4)`

For more information please see [reference](https://manual.gromacs.org/current/reference-manual/file-formats.html#gro).

---------------------------------------------------------------------

## `gro_open`
Open file for reading and checked. Frame offset table is constructed (byte position where frame starts). Move to selected frame with [fseek](https://gcc.gnu.org/onlinedocs/gfortran/FSEEK.html#FSEEK). If the function fails stat is a positive number, otherwise it is set to zero. If stat is not present end error occurs error stop is raised.
```fortran
call gro_open (file, unit, nframes, offset, stat, errmsg)
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

## `gro_close`
Close opened file.
```fortran
call gro_close (unit, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `glo_read_header`
Read file header in GRO format.
```fortran
call gro_read_header (unit, comment, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`character(*), intent(OUT) :: comment`
> Comment string. Optional time in ps after `t=`.

`integer, intent(OUT) :: natoms`
> Number of atoms.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `gro_read_data`
Read data in GRO format.
```fortran
call gro_read_data (unit, natoms, atomi, atomn, resi, resn, coor, box, vel, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, intent(OUT) :: atomi, resi`
> Atom and residue index.

`character(*), intent(OUT) :: atomn, resn`
> Atom and residue name.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in nm.

`real, dimension(DIM,DIM), intent(OUT) :: box`
> Simulation box dimensions in nm.

`real, dimension(DIM,natoms), intent(OUT), optional :: vel`
> Atom velocity in nm/ps.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `gro_read_coor`
Read ONLY coordinate data in GRO format.
```fortran
call gro_read_coor (unit, natoms, coor, box, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in nm.

`real, dimension(DIM,DIM), intent(OUT) :: box`
> Simulation box dimensions in nm.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `gro_skip_data`
Read through file (skip) in GRO format but do not store any data.
```fortran
call gro_skip_data (unit, natoms, stat, errmsg)
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

## `gro_write`
Write data in GRO format.
```fortran
call gro_write (unit, comment, natoms, atomi, atomn, resi, resn, coor, box, vel, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`character(*), intent(OUT) :: comment`
> Comment string.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, intent(IN) :: atomi, resi`
> Atom and residue index.

`character(*), intent(IN) :: atomn, resn`
> Atom and residue name.

`real, dimension(DIM,natoms), intent(IN) :: coor`
> Atom cartesian coordinates in nm.

`real, dimension(DIM,DIM), intent(IN) :: box`
> Simulation box dimensions in nm.

`real, dimension(DIM,natoms), intent(IN), optional :: vel`
> Atom velocity in nm/ps.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## Example
Example of reading GRO file.
```fortran
program main
  use atomlib_groio
  use iso_fortran_env, only: INT64
  implicit none
  integer, parameter :: DIM = 3
  integer, allocatable :: atomi(:), resi(:)
  character(6), allocatable :: atomn(:), resn(:)
  real, allocatable :: coor(:,:)
  integer(INT64), allocatable :: offset(:)
  character(128) :: comment
  real :: box(DIM,DIM)
  integer :: unit, natoms, nframes, status

  call gro_open ("conf.gro", unit, nframes, offset, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call gro_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  allocate (atomi(natoms), resi(natoms), atomn(natoms), &
  & resn(natoms), coor(DIM,natoms), STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call gro_read_data (unit, natoms, atomi, atomn, resi, resn, coor, box, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call gro_close (unit, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message  

end program main
```