####################################################################################################
### A8_computation.g
### 
### This program computes the order of the solution group Γ(M_A8, 0), which shows that the solution
### group is isomorphic to A8.
### 
### The estimated total runtime of this program is:  between 50 and 60 hours (!).
### The estimated memory usage of this program is:   between 50 and 60 GB RAM.
### 
####################################################################################################
# Copyright (C) 2024  J. van Dobben de Bruyn, Technical University of Denmark.
# 
# SPDX-License-Identifier: GPL-2.0-or-later
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
####################################################################################################



Read(Filename(DirectoryCurrent(),"build_system_from_group.g"));
Read(Filename(DirectoryCurrent(),"solution_group.g"));

Print("Testing the alternating group A8...\n");

# Compute the matrix M_A8.
A8 := AlternatingGroup(8);
M_A8 := build_system_from_group(A8);

# Check that the matrix has the right dimensions and rank (over GF(2)).
Assert(0, DimensionsMat(M_A8) = [1645, 315]);
Assert(0, RankMat(Identity(GF(2)) * M_A8) = 315);

# Compute the homogeneous solution group.
b := 0 * [1..1645];
G := solution_group(M_A8, b);

# Check that the solution group is perfect and compute its order.
Assert(0, IsPerfect(G));
Print("Computing the order of the solution group Gamma(M_A8). This will take around 50 to 60 hours...\n");
o := Order(G);
Print("Gamma(M_A8) is a perfect group of order ", o, ".\n");
Assert(0, o = 20160);

# This identifies Gamma(M_A8) uniquely, since A8 is the only group
# of order |G| = |A8| = 20160 that has A8 as a quotient.


