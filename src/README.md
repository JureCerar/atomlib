# TO-DO:
- [ ] As reference use PyMol exposed types (some of them): https://pymolwiki.org/index.php/Iterate
- [ ] Combine all structure and trajcetory types in single type.
- [ ] Implement NAMD .coor data type: https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2009-2010/4104.html
- [ ] Add error number and detailed error string.
- [ ] Add CHEM and ATOM modules with fuctions.
- [ ] (?) List through Python math & string library for ideas.
- [ ] ...
- [x] Added sorting algorithm `xslib_sort` (Not part of main library).
- [ ] Fix atom/res in gro files.
- [ ] Add `xslib_algorithms` functions (qsort, nearest, linest, minimize, etc.).
- [ ] Add `xslib_plot` to main library (+rework)
- [ ] Add create mask(resn, resi, etc) (tpl->gro) for `tpl_t`.
- [ ] Rework csv_t file
- [ ] Add allocatable atribute to fileio. 
- [ ] Add constructors to data files
- [ ] Add operations(+) to coordinate files.
- [ ] Add copy I/O peration to coordinate files.
- [ ] Add initialize functions to data files.
* obj = obj + obj_t( xyz(:), "name" )

## Notes
- Trajectory read also structure file.
'''
stat = traj%read( file,  first, last, stride, structure=str )
'''