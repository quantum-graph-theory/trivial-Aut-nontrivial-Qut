####################################################################################################
### Mathieu.g
### 
### This program builds the linear system M_H (defined in [DSR23+, Definition 6.1]) where H is one
### of the Mathieu groups M9, M10, M11, M12, M21, M22, M23, M24 and checks whether or not this
### system has full column rank.
### 
### The system M_H has full column rank if and only if the homogeneous solution group Γ(M_H, 0) is
### perfect, so this program tests for which Mathieu groups H the homogeneous solution group
### Γ(M_H, 0) is perfect. For brevity, we say that such groups H "work".
### 
### This program finds the following results:
### 
###     * The five sporadic Mathieu groups M11, M12, M22, M23 and M24 work;
###     
###     * The other, non-sporadic Mathieu groups M9, M10 ≤ M11 and M21 ≤ M22 do not work.
### 
### The estimated total runtime of this program is:  between 45 and 60 minutes.
### The estimated memory usage of this program is:   between 60 and 70 GB RAM.
### 
### 
### In addition to the groups found by this program, we have the following results:
### 
###     * The proof in our paper [DRS23+] shows that the alternating group A_n works for all n ≥ 7.
###     
###     * The program "small_perfect_groups.g" finds all perfect groups of order ≤ 10080 that work.
###     
###     * The program "PSL2q.g" shows that the projective special linear group PSL(2, q) works for
###       all odd prime powers 19 ≤ q ≤ 349. This suggests that PSL(2, q) works for all odd prime
###       powers q ≥ 19, though we don't have a proof.
### 
### 
### References:
### 
###       [DRS23+].  Josse van Dobben de Bruyn, David E. Roberson, and Simon Schmidt.
###                  “Asymmetric graphs with quantum symmetry”.
###                  Preprint, 2023.  https://arxiv.org/abs/2311.04889
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



Read(Filename(DirectoryCurrent(),"build_sparse_F2_system_from_group.g"));

# Make and print the list of the "degrees" of the Mathieu groups available in GAP;
# see https://docs.gap-system.org/doc/ref/chap50.html#X788FA7DE84E0FE6A
idx_list := [9, 10, 11, 12, 21, 22, 23, 24];
Print("Testing the Mathieu groups");
for k in idx_list do
	if k = 24 then
		Print(", and");
	elif k > 9 then
		Print(",");
	fi;
	Print(" M_{", k, "}");
od;
Print("...\n\n");
successful_ids := [];

for k in idx_list do
	Mk := MathieuGroup(k);
	Print("\nBuilding the linear system M_{M", k, "} for the Mathieu group M", k, " (as a sparse matrix)");
	if k = 24 then
		Print(". This will take around 10 to 15 minutes");
	fi;
	Print("...\n");
	M_Mk := build_sparse_F2_system_from_group(Mk);
	num_vars := Ncols(M_Mk);
	num_eqs := Nrows(M_Mk);
	Print("Found ", num_vars , " variables and ", num_eqs, " equations.\n");
	
	Print("Computing the rank of the system M_{M", k, "}");
	if k = 24 then
		Print(". This will take around 25 to 35 minutes");
	fi;
	Print("...\n");
	r := Rank(TransposedSparseMat(M_Mk)); # with these linear systems, computing the rank turns out to be MUCH faster after first taking the transpose
	Print("Found rank(M_{M", k, "}) = ", r, ".\n");
	if r = num_vars then
		Print("Conclusion: M_{M", k, "} has full column rank, so the solution group Γ(M_{M", k, "}, 0) is a non-trivial perfect group. Verdict: Success!\n");
		Append(successful_ids, [k]);
	else
		Print("Conclusion: M_{M", k, "} does not have full column rank, so the solution group Γ(M_{M", k, "}, 0) is not a perfect group. Verdict: fail.\n");
	fi;
od;

Print("\nSummary: the following Mathieu groups give rise to a linear system M_H");
Print("\n         with full column rank, so that the solution group Γ(M_H, 0)");
Print("\n         is a non-trivial perfect group:\n");
for id in successful_ids do
	Print("    M_{", id, "}\n");
od;


