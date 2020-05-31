#========================================================================
# This file is part of EwE-F
# Copyright (C) 2011-2019 METU Institute of Marine Sciences, Turkey.
#
# This program is free software; you can redistribute it and/or modify 
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with  this program; if not, 
# see <http://www.gnu.org/licenses/gpl-2.0.html>.
#========================================================================
# begin

# Compiler definitions and parameters

# Set compiler:  Intel Fortran Compiler (ifort)
# Oracle Solaris Studio Fortran Compiler (g95) or
# GNU Fortran Compiler (gfortran) 
FC = gfortran

# Uncomment below to build for 64 bits architecture
# FCFLAGS = -m64

# Uncomment below for double precison
# CPPDEFS = -D_isDP_

# Set the path to your HDF5 installation
# The HDF5 should have been compiled by the same compiler
# specified above
HDF5DIR = /usr/local/hdf5-1.10.4

# This flag below specifies the directory into which the compiled module
# files (*.mod) are going to be put. For GNU Fortran compiler
# set it to "-J", for Oracle Solaris Studio Fortran Compiler "-moddir=" and
# for Intel Fortran and PGI Fortran compilers set it to  "-module " (with a whitespace afterwards!!)
MODFLAG = -J

# EwEFORTRAN root directory
EwEFORTRANDIR = $(CURDIR)

# Directory names for the source, object, module and binary files
SRCDIR = $(EwEFORTRANDIR)/src
OBJDIR = $(EwEFORTRANDIR)/obj
MODDIR = $(EwEFORTRANDIR)/mod
BINDIR = $(EwEFORTRANDIR)/bin

# Compiler and linker flags for Intel Fortran compiler
# CFLAGS = -c -g -CB -traceback -warn all -warn errors $(FCFLAGS)
# LFLAGS = -g -CB -traceback -warn all -warn errors $(FCFLAGS)

# Compiler and linker flags for PGI Fortran compiler
# CFLAGS = -c -g -Mbounds -traceback -w $(FCFLAGS)
# LFLAGS = -g -Mbounds -traceback -w $(FCFLAGS)

# Compiler and linker flags for GNU Fortran compiler
 CFLAGS = -c -g -Warray-bounds -fbacktrace -fbounds-check -Wall $(FCFLAGS)
 LFLAGS = -g -Warray-bounds -fbacktrace -fbounds-check -Wall $(FCFLAGS)

# Compiler and linker flags for Oracle Solaris Studio Fortran compiler
# CFLAGS = -c -g -w3 -traceback
# LFLAGS = -g -w3 -traceback

####### The section below this notice does not need any changes ########

# Flags for module directory
CFLAGS += $(MODFLAG)$(MODDIR)
LFLAGS += $(MODFLAG)$(MODDIR)

# HDF5 flags
CFLAGS += -I$(HDF5DIR)/include -L$(HDF5DIR)/lib -lhdf5_fortran -lhdf5
LFLAGS += -I$(HDF5DIR)/include -L$(HDF5DIR)/lib -lhdf5_fortran -lhdf5

# Source files
SRCECOPATH    = $(SRCDIR)/ecopath.F90
SRCMODECOPATH = $(SRCDIR)/statevartypesecopath.F90
SRCMODHDF5    = $(SRCDIR)/readHDF5Database.F90
SRCGROUP      = $(SRCDIR)/readEcopathScenario_io.F90
SRCDET        = $(SRCDIR)/readDetritusFate_io.F90
SRCDIET       = $(SRCDIR)/readDietComposition_io.F90
SRCORDER      = $(SRCDIR)/calculateSolvingOrderOfLinearEquations.F90
SRCENERGY     = $(SRCDIR)/calculateEnergyBalance.F90
SRCSANITY     = $(SRCDIR)/warnSanityChecks.F90
SRCHDF5       = $(SRCDIR)/writeHDF5ResultsFile.F90
SRCALLOC      = $(SRCDIR)/setEcopathScenarioParameters.F90
SRCCONSUME    = $(SRCDIR)/calculatePredation.F90
SRCFLOWS2DET  = $(SRCDIR)/calculateDetritalFlows.F90
SRCCALCDF     = $(SRCDIR)/calculateDetritusFate.F90
SRCCALCBADET  = $(SRCDIR)/calculateBAofDetritus.F90
SRCCALCEEDET  = $(SRCDIR)/calculateEEofDetritus.F90
SRCSTANZA     = $(SRCDIR)/calculateMultistanzaParameters.F90

# Object files
OBJECOPATH    = $(OBJDIR)/ecopath.o
OBJMODECOPATH = $(OBJDIR)/statevartypesecopath.o
OBJMODHDF5    = $(OBJDIR)/readHDF5Database.o
OBJGROUP      = $(OBJDIR)/readEcopathScenario_io.o
OBJDET        = $(OBJDIR)/readDetritusFate_io.o
OBJDIET       = $(OBJDIR)/readDietComposition_io.o
OBJORDER      = $(OBJDIR)/calculateSolvingOrderOfLinearEquations.o
OBJENERGY     = $(OBJDIR)/calculateEnergyBalance.o
OBJSANITY     = $(OBJDIR)/warnSanityChecks.o
OBJHDF5       = $(OBJDIR)/writeHDF5ResultsFile.o
OBJALLOC      = $(OBJDIR)/setEcopathScenarioParameters.o
OBJCONSUME    = $(OBJDIR)/calculatePredation.o
OBJFLOWS2DET  = $(OBJDIR)/calculateDetritalFlows.o
OBJCALCDF     = $(OBJDIR)/calculateDetritusFate.o
OBJCALCBADET  = $(OBJDIR)/calculateBAofDetritus.o
OBJCALCEEDET  = $(OBJDIR)/calculateEEofDetritus.o
OBJSTANZA     = $(OBJDIR)/calculateMultistanzaParameters.o

OBJECTS = $(OBJMODECOPATH) $(OBJMODHDF5) $(OBJHDF5) $(OBJALLOC) \
$(OBJCONSUME) $(OBJFLOWS2DET) $(OBJCALCDF) $(OBJCALCBADET) $(OBJCALCEEDET) $(OBJSTANZA) \
$(OBJGROUP) $(OBJDET) $(OBJDIET) $(OBJORDER) $(OBJENERGY) $(OBJECOPATH) \
$(OBJSANITY)

# Module files
MODECOPATH = $(MODDIR)/statevartypesecopath.mod
MODHDF5 = $(MODDIR)/readhdf5database.mod

MODULES = $(MODECOPATH) $(MODHDF5)

# Binary files
BINECOPATH     = $(BINDIR)/ecopath.out

# Now compile the program
$(BINECOPATH): $(OBJECTS)
	$(FC) $(CPPDEFS) $(LFLAGS) $(OBJECTS) -o $(BINECOPATH)

$(OBJMODECOPATH): $(SRCMODECOPATH)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCMODECOPATH) -o $(OBJMODECOPATH)

$(OBJMODHDF5): $(SRCMODHDF5)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCMODHDF5) -o $(OBJMODHDF5)

$(OBJGROUP): $(SRCGROUP)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCGROUP) -o $(OBJGROUP)

$(OBJDET): $(SRCDET)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCDET) -o $(OBJDET)

$(OBJDIET): $(SRCDIET)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCDIET) -o $(OBJDIET)

$(OBJORDER): $(SRCORDER)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCORDER) -o $(OBJORDER)

$(OBJENERGY): $(SRCENERGY)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCENERGY) -o $(OBJENERGY)

$(OBJSANITY): $(SRCSANITY)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCSANITY) -o $(OBJSANITY)

$(OBJHDF5): $(SRCHDF5)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCHDF5) -o $(OBJHDF5)

$(OBJALLOC): $(SRCALLOC)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCALLOC) -o $(OBJALLOC)

$(OBJCONSUME): $(SRCCONSUME)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCCONSUME) -o $(OBJCONSUME)

$(OBJFLOWS2DET): $(SRCFLOWS2DET)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCFLOWS2DET) -o $(OBJFLOWS2DET)

$(OBJCALCDF): $(SRCCALCDF)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCCALCDF) -o $(OBJCALCDF)

$(OBJCALCBADET): $(SRCCALCBADET)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCCALCBADET) -o $(OBJCALCBADET)

$(OBJCALCEEDET): $(SRCCALCEEDET)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCCALCEEDET) -o $(OBJCALCEEDET)

$(OBJSTANZA): $(SRCSTANZA)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCSTANZA) -o $(OBJSTANZA)

$(OBJECOPATH): $(MODULES) $(SRCECOPATH)
	$(FC) $(CPPDEFS) $(CFLAGS) $(SRCECOPATH) -o $(OBJECOPATH)

clean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(BINECOPATH)

# end
