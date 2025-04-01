####################################################################################################
### build_sparse_F2_system_from_group.g
### 
### This file provides a function that, given a group H, constructs the matrix M_H that encodes all
### triples of pairwise commuting elements of order 2 in H that multiply to 1.
### 
### The difference between this file and the file "build_system_from_group.g" is that this file
### returns the system as a sparse matrix over F2, using a data type from the GAP package Gauss,
### whereas the file "build_system_from_group.g" returns an ordinary GAP matrix with integer
### entries, which has a much higher memory footprint (especially for large groups).
### 
### The sparse version provided in this file is necessary to verify that our construction can also
### be carried out starting from the Mathieu group M24 and other large groups, but for small groups
### the version from "build_system_from_group.g" suffices.
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



LoadPackage("Gauss");


build_sparse_F2_system_from_group := function(G)
	local order_2_elms, l, index_of_element, rels, r, i, j, p, k, cur_eq, M;
	
	# Find all elements of order 2 in G.
	order_2_elms := Filtered(Elements(G), x -> Order(x) = 2);
	l := Size(order_2_elms);
	
	# Part of the remaining code assumes that there is at least one element of order 2,
	# so we catch the corner case of G not having any elements of order 2.
	if l = 0 then
		return [];
	fi;
	
	# Create a dictionary that stores the index of every element of order 2 in the list.
	index_of_element := NewDictionary(order_2_elms[1], true, order_2_elms);
	for i in [1 .. l] do
		Assert(0, LookupDictionary(index_of_element, order_2_elms[i]) = fail);
		AddDictionary(index_of_element, order_2_elms[i], i);
		Assert(0, LookupDictionary(index_of_element, order_2_elms[i]) = i);
	od;
	
	# Find all relations between three pairwise commuting element of order 2.
	# We do this by generating all pairs of elements of order 2. If they commute, then their
	# product is also an element of order 2. We look it up in the dictionary, compare its index k
	# to the indices i and j of the generated pair (to avoid duplicate entries), and if we have
	# i < j < k then we add an equation corresponding to this triple to the linear system.
	rels := [];
	r := 0;
	for i in [1 .. l] do
		for j in [i + 1 .. l] do
			if Comm(order_2_elms[i], order_2_elms[j]) = Identity(G) then
				p := order_2_elms[i] * order_2_elms[j];
				Assert(0, p = order_2_elms[j] * order_2_elms[i]);
				k := LookupDictionary(index_of_element, p);
				Assert(0, k <> fail);
				Assert(0, Comm(order_2_elms[i], order_2_elms[k]) = Identity(G));
				Assert(0, Comm(order_2_elms[j], order_2_elms[k]) = Identity(G));
				Assert(0, order_2_elms[i] * order_2_elms[j] * order_2_elms[k] = Identity(G));
				if k > j then
					r := r + 1;
					Append(rels, [[i, j, k]]);
				fi;
			fi;
		od;
	od;
	# Build the matrix from the list of lists of indices.
	# Note: since we only specify the indices of the non-zero entries of the matrix without specifying
	# the corresponding values, the function SparseMatrix will infer that the ground field is F2.
	# Therefore the returned sparse matrix will be a matrix over F2.
	# (For more information, see the documentation for the GAP package Gauss.)
	Assert(0, Size(rels) = r);
	M := SparseMatrix(r, l, rels);;
	
	# Return the matrix we found.
	return M;
	end;;


