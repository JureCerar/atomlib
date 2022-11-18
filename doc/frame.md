# frame

## `frame_t`
Trajectory frame.
```fortran
  type :: frame_t
    real, allocatable :: coor(:,:)
    real :: box(DIM,DIM)
    real :: time = 0.0
    integer :: natoms = 0
  contains
    procedure :: allocate 
    procedure :: reallocate 
    procedure :: deallocate 
    procedure :: allocated
    generic :: operator(+)
    generic :: assignment(=)
    generic :: write(formatted)
  end type frame_t  
```
#### *Values*
`real, allocatable, dimension(:,:) :: coor`
> Atom's cartesian coordinates (x,y,z) in nm.

`real, dimension(DIM,DIM) :: box`
> Simulation box size vector in nm.

`real :: time`
> Frame time in ps.

`integer :: natoms`
> Number of atoms in configuration.

---------------------------------------------------------------------

## `constructor`
Constructor for `frame_t` derived type.
```
out = frame_t(natoms, time, coor, box)
```
#### *Parameters*

---------------------------------------------------------------------

## `allocate`
Allocate object to specified number of atoms.
```fortran
call frame%allocate (natoms, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `reallocate`
Reallocate object to specified number of atoms (bigger or smaller).
```fortran
call frame%reallocate (natoms, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `deallocate`
Deallocate object.
```fortran
call frame%deallocate (stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `allocated`
Check if object is (properly) allocated.
```fortran
out = frame%allocated()
```
#### *Parameters*

---------------------------------------------------------------------

## `assignment(=)`
Assign value to object. Delete data if already allocated.
```fortran
out = frame
```
#### *Parameters*

---------------------------------------------------------------------

## `operator(+)`
Append data to another object. Both object must be initialized.
```fortran
out = frame + frame
```
#### *Parameters*

---------------------------------------------------------------------

## `write(formatted)`
Formatted write to output unit. Use `DT` or `DT(w,d)`use for custom output format where `w` is the number of positions to be used and `d` is the number of digits to the right of the decimal point.
```fortran
write (*,*) frame
```
