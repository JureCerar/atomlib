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

set ( EXPORT_INSTALL_DIR "${INSTALL_LIBDIR}/cmake/${PACKAGE_VERSION}" )

include ( CMakePackageConfigHelpers ) # Standard CMake module
write_basic_package_version_file (
  "${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake"
  VERSION ${VERSION}
  COMPATIBILITY SameMajorVersion
)

# Install package config file
configure_package_config_file (
  "${CMAKE_SOURCE_DIR}/cmake/pkg/${CMAKE_PROJECT_NAME}-config.cmake.in"
  "${CMAKE_BINARY_DIR}/pkg/${PACKAGE_NAME}-config.cmake"
  INSTALL_DESTINATION "${EXPORT_INSTALL_DIR}"
  PATH_VARS EXPORT_INSTALL_DIR
)

# Install the config and version files so that we can find this project with others
install ( FILES
  "${CMAKE_BINARY_DIR}/pkg/${PACKAGE_NAME}-config.cmake"
  "${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake"
  DESTINATION "${EXPORT_INSTALL_DIR}"
)
configure_file(
 "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${CMAKE_PROJECT_NAME}.pc.cmake.in"
 "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}.pc"
 @ONLY
)
install ( FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}.pc"
  DESTINATION "${LIBDIR}/pkgconfig"
)
