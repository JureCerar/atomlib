# trajio

## `traj_t`
```fortran
  type :: traj_t
    type(frame_t), allocatable :: frame(:)
    integer :: nframes
  contains
    procedure :: allocate
    procedure :: reallocate
    procedure :: deallocate
    procedure :: allocated
    generic :: write(formatted)
    procedure :: open
    procedure :: close
    procedure :: fseek
    procedure :: ftell
    procedure :: read_next
    procedure :: save
    procedure :: load
  end type traj_t
```
#### *Values*
`type(frame_t), allocatable, dimension(:) :: frame`
> Configuration frame.

`integer :: natoms`
> Number of frames in trajectory.

---------------------------------------------------------------------

## `constructor`
Constructor for `traj_t` derived type.
```
out = traj_t(nframes, frame)
```
#### *Parameters*

---------------------------------------------------------------------

## `allocate`
Allocate object to specified number of frames.
```fortran
call obj%allocate (nframes, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `reallocate`
Reallocate object to specified number of frames (bigger or smaller).
```fortran
call obj%reallocate (nframes, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `deallocate`
Deallocate object.
```fortran
call obj%deallocate (stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `allocated`
Check if object is (properly) allocated.
```fortran
out = obj%allocated()
```
#### *Parameters*

---------------------------------------------------------------------

## `write(formatted)`
Formatted write to output unit. Use `DT` or `DT(w,d)`use for custom output format where `w` is the number of positions to be used and `d` is the number of digits to the right of the decimal point.
```fortran
write (*, *) obj
```
#### *Parameters*

---------------------------------------------------------------------

## `open`
Open trajectory for reading.
```fortran
call obj%open (file, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `close`
Close trajectory.
```fortran
call obj%close (stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `fseek`
Move position to the specified frame from whence.
```fortran
call obj%fseek (offset, whence, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `ftell`
Returns current position in trajectory file.
```fortran
out = obj%ftell()
```
#### *Parameters*

---------------------------------------------------------------------

## `read_next`
Read next frame (or next N frames) in trajectory file. 
```fortran
call out%read_next (nframes, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `load`
Load (read) trajectory file. File format is determined from file extension.
```fortran
call obj%load (file, first, last, stride, stat, errmsg)
```
#### *Parameters*

---------------------------------------------------------------------

## `save`
Save (write) configuration file. File format is determined from file extension.
```fortran
call obj%save (file, stat, errmsg)
```
#### *Parameters*
