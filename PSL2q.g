####################################################################################################
### PSL2q.g
### 
### This program builds the linear system M_H (defined in [DSR23+, Definition 6.1]), where H is a
### projective special linear group PSL(2, q) for some prime power q, and checks whether or not this
### system has full column rank.
### 
### The system M_H has full column rank if and only if the homogeneous solution group Γ(M_H, 0) is
### perfect, so this program finds all groups of the form H = PSL(2, q) for some prime power q
### (up to a certain maximum) for which the homogeneous solution group Γ(M_H, 0) is perfect.
### For brevity, we say that such groups H "work".
### 
### This program finds the following results:
### 
###     * The projective special linear group PSL(2, q) works whenever q is odd and 19 ≤ q ≤ 349;
###     
###     * The projective special linear group PSL(2, q) does not work if q is even or q ≤ 17.
### 
### We did not test for q > 350. We suspect that PSL(2, q) works whenever q is odd and q ≥ 19,
### though we don't have a proof.
### 
### The estimated total runtime of this program is:  between 45 and 50 hours (!).
### The estimated memory usage of this program is:   between 100 and 120 GB RAM.
### 
### 
### In addition to the groups found by this program, we have the following results:
### 
###     * The proof in our paper [DRS23+] shows that the alternating group A_n works for all n ≥ 7.
###     
###     * The program "small_perfect_groups.g" finds all perfect groups of order ≤ 10080 that work.
###     
###     * The program "Mathieu.g" shows that the sporadic Mathieu groups M11, M12, M22, M23 and M24
###       also work, whereas the non-sporadic Mathieu subgroups M9, M10 ≤ M11 and M21 ≤ M22 do not.
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

max_q := 350;

Print("Testing projective special linear groups PSL(2, q) for all prime powers q ≤ ", max_q, "...\n\n");
successful_ids := [];

for q in [1..max_q] do
	if IsPrimePowerInt(q) then
		H := PSL(2, q);
		Print("Trying PSL(2, ", q, ")... \c");
		M_H := build_sparse_F2_system_from_group(H);
		if M_H = [] then
			Print("Fail.\n");
			continue;
		fi;
		num_vars := Ncols(M_H);
		num_eqs := Nrows(M_H);
		r := Rank(TransposedSparseMat(M_H)); # with these linear systems, computing the rank turns out to be MUCH faster after first taking the transpose
		if r = num_vars then
			Print("Success!\n");
			Append(successful_ids, [q]);
		else
			Print("Fail.\n");
		fi;
	fi;
od;

Print("\nSummary: the following projective special linear groups PSL(2, q) with q ≤ ", max_q);
Print("\n         give rise to a linear system M_H with full column rank, so that");
Print("\n         the solution group Γ(M_H, 0) is a non-trivial perfect group:\n");
for q in successful_ids do
	Print("    PSL(2, ", q, ")\n");
od;


