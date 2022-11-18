# This file is part of atomlib
# https://github.com/JureCerar/atomlib
#
# Copyright (C) 2022 Jure Cerar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Include library and module files from the build
link_directories ( "${CMAKE_BINARY_DIR}/lib" )
include_directories ( "${CMAKE_BINARY_DIR}/include" )

# Get test source files and copy sample files.
file ( GLOB "TEST_SRC" RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.f90" "*.F90" )
file ( COPY "${CMAKE_CURRENT_SOURCE_DIR}/files/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/files/" )

# Make custom `check` compile target.
add_custom_target ( check COMMAND ${CMAKE_CTEST_COMMAND} )
add_dependencies ( check ${CMAKE_PROJECT_NAME}_static )

# Compile tests as part of 'check'.
foreach ( "FILE" ${TEST_SRC} )
  get_filename_component ( "NAME" ${FILE} NAME_WE )
  add_executable ( ${NAME} EXCLUDE_FROM_ALL ${FILE} )
  target_link_libraries ( ${NAME} ${CMAKE_PROJECT_NAME}_static )
  add_dependencies ( check ${NAME} )
endforeach ()

# Find all test files
file ( GLOB "CONF_FILES" "${CMAKE_CURRENT_BINARY_DIR}/files/conf.*" )
file ( GLOB "TRAJ_FILES" "${CMAKE_CURRENT_BINARY_DIR}/files/traj.*" )

# Add tests
add_test ( NAME "atomlib" COMMAND atomlib )
add_test ( NAME "xyzio" COMMAND xyzio "files/conf.xyz" )
add_test ( NAME "groio" COMMAND groio "files/conf.gro" )
add_test ( NAME "pdbio" COMMAND pdbio "files/conf.pdb" )
add_test ( NAME "crdio" COMMAND crdio "files/conf.crd" "files/conf.cor" )
add_test ( NAME "trrio" COMMAND trrio "files/traj.trr" )
add_test ( NAME "xtcio" COMMAND xtcio "files/traj.xtc" )
add_test ( NAME "dcdio" COMMAND dcdio "files/traj.dcd" )
add_test ( NAME "confio" COMMAND confio ${CONF_FILES} )
add_test ( NAME "trajio" COMMAND trajio ${CONF_FILES} ${TRAJ_FILES} )
add_test ( NAME "ndxio" COMMAND ndxio "files/index.ndx" )
add_test ( NAME "pdhio" COMMAND pdhio "files/data.pdh" )
