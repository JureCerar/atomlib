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

# Compilers & linker
FC ?= gfortran
FC_FLAGS := -g -Og -march=native -Wall -fcheck=all -fPIC
CC ?= gcc
CC_FLAGS := -g -Og -march=native -fPIC
# Install
INSTALL := install
INSTALL_DIR := /usr/local
# Source
SRC := xdrfile.c xdrfor.f90 xyzio.f90 groio.f90 pdbio.f90 crdio.f90
SRC += trrio.f90 xtcio.f90 dcdio.f90 frame.f90 confio.f90 trajio.f90 
SRC += ndxio.f90 pdhio.f90 atomlib.f90
# Objects
OBJ := $(SRC)
OBJ := $(OBJ:.F90=.o)
OBJ := $(OBJ:.f90=.o)
OBJ := $(OBJ:.c=.o)
# Library (& binary)
LIB = libatomlib.so
BIN = a.out

# Build
default: $(LIB)

$(LIB): $(OBJ)
	$(FC) -shared $^ -o $(LIB) $(FC_FLAGS)

$(BIN): $(OBJ)
	$(FC) $^ -o $(BIN) $(FC_FLAGS)

%.o: %.f90
	$(FC) -c $< $(FC_FLAGS)

%.o: %.F90
	$(FC) -c $< $(FC_FLAGS)

%.o: %.c
	$(CC) -c $< $(CC_FLAGS)

.PHONY: test
test: $(BIN)
	@./$(BIN)

.PHONY: install
install:
	$(INSTALL) $(LIB) $(INSTALL_DIR)/lib
	$(INSTALL) *.mod  $(INSTALL_DIR)/include

.PHONY: clean
clean:
	@rm -rf *.o *.mod $(BIN) $(LIB)

.PHONY: help	
help:
	@$(foreach V,$(sort $(.VARIABLES)), \
	$(if $(filter-out environment% default automatic, \
	$(origin $V)),$(warning $V=$($V) ($(value $V)))))