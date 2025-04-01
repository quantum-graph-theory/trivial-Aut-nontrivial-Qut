####################################################################################################
### A7_computation.g
### 
### This program computes the order of the solution group Γ(M_A7, 0), which shows that the solution
### group is isomorphic to 3.A7 (the triple cover of A7).
### 
### The estimated total runtime of this program is:  between 30 and 60 minutes.
### The estimated memory usage of this program is:   less than 4 GB RAM.
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

Print("Testing the alternating group A7...\n");

# Compute the matrix M_A7.
A7 := AlternatingGroup(7);
M_A7 := build_system_from_group(A7);

# Check that the matrix has the right dimensions and rank (over GF(2)).
Assert(0, DimensionsMat(M_A7) = [140, 105]);
Assert(0, RankMat(Identity(GF(2)) * M_A7) = 105);

# Compute the homogeneous solution group.
b := 0 * [1..140];
G := solution_group(M_A7, b);

# Check that the solution group is perfect and compute its order.
Assert(0, IsPerfect(G));
Print("Computing the order of the solution group Gamma(M_A7). This will take around 30 to 60 minutes...\n");
o := Order(G);
Print("Gamma(M_A7) is a perfect group of order ", o, ".\n");
Assert(0, o = 7560);

# This identifies Gamma(M_A7) uniquely, since there is exactly one perfect
# group of order 7560; namely, the group 3.A7 (the triple cover of A7).


