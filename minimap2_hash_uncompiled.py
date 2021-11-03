#!/usr/bin/env python3
"""
Wrapper for minimap2's hash function.

This is an uncompiled version to use if minimap2_hash.pyx is unavailable or
can't be compiled.

The hash function, minimap2_hash(), is derived from the minimap2 source code
(the function hash64 in sketch.c). It is modified here by the incorporation of
a seed. If the seed is zero, the hash is identical to the minimap2 function. 

For a given seed and a mask of the form (4**k)-1, minimap2_hash(seed,*,mask) is
invertible over 4**k.
"""

# minimap2_hash--
#	seed: a 64-bit seed for the hash function
#	v:    the 64-bit value to hash
#	mask: the bits to hash; usually this is (4**k)-1 where k is the kmer size

def minimap2_hash(seed,v,mask):  # (all arguments are 64-bit values)
	u = (v + seed) & mask
	u = ((~u) + (u << 21)) & mask			# u = (u<<21)-(u+1) = 77594587*u-1
	u = u ^ (u >> 24)
	u = ((u + (u << 3)) + (u << 8)) & mask  # u *= 265
	u = u ^ (u >> 14)
	u = (u + (seed >> 5)) & mask
	u = ((u + (u << 2)) + (u << 4)) & mask  # u *= 21
	u = u ^ (u >> 28)
	u = (u + (u << 31)) & mask				# u *= 2147483649
	return u
