# confio
<!-- One object to rule them all, one object to find them,
  One object to bring them all, and in the darkness bind them -->

## `conf_t`
`conf_t` is object designed to handle all molecular cinfiguration files.
```fortran
type :: conf_t
  integer, allocatable           :: atomi(:)
  character(MAXLEN), allocatable :: atomn(:)
  integer, allocatable           :: resi(:)
  character(MAXLEN), allocatable :: resn(:)
  character(MAXLEN), allocatable :: chain(:)
  character(MAXLEN), allocatable :: element(:)
  real, allocatable              :: bfact(:)
  real, allocatable              :: pcharge(:)
  integer, allocatable           :: charge(:)
  real, allocatable              :: coor(:,:)
  real                           :: box(DIM,DIM) = 0.0
  integer                        :: natoms = 0
contains
  procedure :: allocate
  procedure :: reallocate
  procedure :: deallocate
  procedure :: allocated
  generic   :: operator(+)
  generic   :: assignment(=)
  generic   :: write(formatted)
  procedure :: load
  procedure :: save
end type conf_t
```
#### *Values*
`integer, allocatable, dimension(:) :: atomi, resi`
> Atom and residue index.

`character(*), allocatable, dimension(:) :: atomn, resn`
> Atom and residue name.

`character(*), allocatable, dimension(:) :: chain`
> Peptide chain name.

`character(*), allocatable, dimension(:) :: element`
> Element name.

`real, allocatable, dimension(:) :: bfact`
> Beta (temperature) factor.

`real, allocatable, dimension(:) :: pcharge`
> Atom's partial charge.

`integer, allocatable, dimension(:) :: charge`
> Atom's formal charge.

`real, allocatable, dimension(:,:) :: coor`
> Atom's cartesian coordinates (x,y,z).

`real, dimension(DIM,DIM) :: box`
> Simulation box size vector.

`integer :: natoms`
> Number of atoms in configuration.

---------------------------------------------------------------------

## `constructor`
Constructor for `conf_t` derived type.
```
out = conf_t(natoms, atomi, atomn, resi, resn, chain, element, &
  &  bfact, pcharge, charge, coor, box)
```
#### *Parameters*

---------------------------------------------------------------------

## `allocate`
Allocate object to specified number of atoms.
```fortran
call conf%allocate (natoms, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `reallocate`
Reallocate object to specified number of atoms (bigger or smaller).
```fortran
call conf%reallocate (natoms, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `deallocate`
Deallocate object.
```fortran
call conf%deallocate (stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `allocated`
Deallocate object.
```fortran
out = conf%allocated()
```
#### *Parameters*

---------------------------------------------------------------------

## `assignment(=)`
Assign value to configuration. Dump data if already allocated.
```fortran
out = conf
```
#### *Parameters*

---------------------------------------------------------------------

## `operator(+)`
Append data to another object. Both object must be initialized.
```fortran
out = conf + conf
```
#### *Parameters*

---------------------------------------------------------------------

## `write(formatted)`
Formatted write to output unit. Use `DT` or `DT(w,d)`use for custom output format where `w` is the number of positions to be used and `d` is the number of digits to the right of the decimal point.
```fortran
write (*, *) conf
```
#### *Parameters*

---------------------------------------------------------------------

## `load`
Load (read) configuration file. File format is determined from file extension.
```fortran
call conf%load (file, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `save`
Save (write) configuration file. File format is determined from file extension.
```fortran
call conf%save (file, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## example
```Fortran
program main
  use atomlib_confio
  implicit none
  type(conf_t) :: conf
  character(128) :: message
  integer :: status

  call conf%load ("input.gro", STAT=status, ERMSG=message)
  if (status /= 0) error stop message

  call conf%save ("output.gro", STAT=status, ERMSG=message)
  if (status /= 0) error stop message

end program main
```