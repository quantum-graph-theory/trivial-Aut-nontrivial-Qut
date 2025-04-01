####################################################################################################
### solution_group.g
### 
### This file provides a function that constructs the solution group associated to a given binary
### linear constraint system Ax = b.
### 
### The code in this file is a modification of the code from:
### 
### [PRSS22]. Connor Paddock, Vincent Russo, Turner Silverthorne, and William Slofstra. Supplementary
###           software for "Arkhipov's theorem, graph minors, and linear system nonlocal games".
###           https://github.com/vprusso/graph_incidence_nonlocal_games
### 
####################################################################################################
# Copyright (C) 2024  C. Paddock, V. Russo, T. Silverthorne, W. Slofstra,
#                     J. van Dobben de Bruyn, Technical University of Denmark.
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





####################################################################################################
## Function auto_detect_homogeneous(b).
## 
## A helper function for the function solution_group, to automatically detect whether or not the
## linear system is homogeneous.
## 
## Input:
## 
##     * b: a vector over GF(2).
## 
## Returns:
## 
##     * True if all entries of b are zero; false otherwise.
## 
####################################################################################################

auto_detect_homogeneous := function(b)
	local non_zero_entries_b, homogeneous, hom_str;
	
	# Check whether or not all entries of b are zero.
	non_zero_entries_b := Filtered(b, x -> x <> 0);
	homogeneous := IsEmpty(non_zero_entries_b);
	
	# Print output message.
	if homogeneous then hom_str := "homogeneous"; else hom_str := "inhomogeneous"; fi;
	Print("Auto-detected that the requested solution group should be ", hom_str, ".\n");
	
	# Return answer.
	return homogeneous;
	end;;




####################################################################################################
## Function solution_group(A, b, [homogeneous]).
## 
## Computes the solution group for the binary linear constraint system Ax = b.
## 
## Inputs:
## 
##     * A: an m x n coefficient matrix with entries either in Z or GF(2);
##       
##     * b: a vector of length m with entries either in Z or GF(2);
##       
##     * homogeneous (optional): a boolean indicating whether or not we should set J = 1.
## 
## Returns:
## 
##     * A finitely presented group encoding the solution group for the binary linear constraint
##       system Ax = b.
## 
## 
## Note: the code for this function was taken from the [PRSS22] file
## 
##       "confluent_rewriting_systems/solution_group_utils/BCS_Presentation.g"
## 
## with the following modifications:
## 
##     * Don't ignore empty rows of the matrix A, but print a warning if the corresponding entry
##       in the vector b is non-zero (because this forces J = 1).
##       
##     * Merge functions "BCS_presentation" and "BCS_NO_j_presentation" into a single function with
##       an optional argument "homogeneous" to avoid code duplication. If this optional argument is
##       omitted, it will be detected automatically by the helper function auto_detect_homogeneous,
##       which tests whether or not all entries of the b vector are zero.
##       
##     * Minor changes in formatting (e.g. error messages, names of variables, whitespace, etc.).
## 
####################################################################################################

solution_group := function(A, b, homogeneous...)
	local hom_str, dims, neqs, nvars, F, J, rels, constraint_elements, prod, gens, i, j, k, l;
	
	# Parse optional argument.
	if Length(homogeneous) > 1 then
		Print("Warning: too many arguments. Ignoring excess arguments.\n");
	fi;
	
	if Length(homogeneous) = 0 then
		# Optional argument omitted; auto-detect.
		homogeneous := auto_detect_homogeneous(b);
	else
		if homogeneous[1] = true or homogeneous[1] = false then
			homogeneous := homogeneous[1];
		else
			Print("Warning: unable to parse optional \"homogeneous\" parameter (read: ", homogeneous, "). Auto-detecting.\n");
			homogeneous := auto_detect_homogeneous(b);
		fi;
	fi;
	
	Assert(0, homogeneous = true or homogeneous = false);
	if homogeneous then hom_str := "homogeneous"; else hom_str := "inhomogeneous"; fi;
	Print("Computing ", hom_str, " solution group presentation.\n");
	
	# Read matrix dimensions.
	dims := DimensionsMat(A);
	neqs := dims[1];  # number of equations (rows of the matrix A)
	nvars := dims[2]; # number of variables (columns of the matrix A)
	
	# Include extra generator for the J element.
	F := FreeGroup(nvars + 1);
	gens := GeneratorsOfGroup(F);
	J := gens[nvars + 1];
	
	# Construct the list of relations.
	# First, make the generators have order 2.
	rels := List(gens, g -> g^2);
	
	# If homogeneous is set to true, force J = 1. Otherwise, make J commute with everything.
	if homogeneous then
		Append(rels, [J]);
	else
		for i in [1 .. nvars] do
			Append(rels, [Comm(gens[i], J)]);
		od;
	fi;

	# Add linear constraints.
	for i in [1 .. neqs] do
		# first find all elements in the constraint
		constraint_elements := [];
		prod := J^-b[i];
		for j in [1 .. nvars] do
			if A[i][j] <> 0 then
				Append(constraint_elements, [gens[j]]);
				prod := prod * gens[j]^A[i][j];
			fi;
		od;
		
		# Add the linear constraint (if it's non-trivial).
		if not IsOne(prod) then
			Append(rels, [prod]);
		fi;
		
		# Add the commuting relations.
		# Some of these might be added more than once, but that's not a big problem because
		# most algorithms in GAP automatically simplify the presentation first.
		for k in [1 .. Length(constraint_elements)] do
			for l in [k + 1 .. Length(constraint_elements)] do
				Append(rels, [Comm(constraint_elements[k], constraint_elements[l])]);
			od;
		od;
		
		# Print a warning message in case the current row of A is zero but the corresponding entry in b is non-zero.
		if constraint_elements = [] and b[i] <> 0 then
			Print("Warning: found zero matrix row with non-zero b entry. This forces J = 1.\n");
		fi;
	od;
	
	# Return the solution group.
	return F / rels;
	end;;



