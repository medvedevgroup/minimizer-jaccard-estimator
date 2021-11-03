#!/usr/bin/env python3
"""
Wrapper for the broccohash hash function.

This is an uncompiled version to use if broccohash.pyx is unavailable or can't
be compiled.

This hash function, broccohash(), is derived from a blog post by David Stafford
[1]. The function described there (which is NOT broccohash) has been used, in
augmented form, in SplitMix64 [2,3]. It is modified here by the incorporation
of a seed. If the seed is zero, the hash is identical to Stafford's function.
However, in that case it hashes zero to zero, which for many uses is not ideal.
A non-zero seed is recommended, and thus the first step of broccohash() xors
the caller's seed with a nothing-up-my-sleeve value, making a non-zero seed
more or less 'automatic'.

Note that Stafford's function (which is called Mix13 in [1]) is similar to the
avalanche step of murmurhash3 [4,5], but allegedly has better mixing properties.

For a given seed, broccohash(seed,*) is invertible over 2**64.

References:
  [1] Stafford, David, "Better Bit Mixing - Improving on MurmurHash3's
      64-bit Finalizer"
      http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html
  [2] Steele Jr, Guy L., Doug Lea, and Christine H. Flood. "Fast splittable
	     pseudorandom number generators." ACM SIGPLAN Notices. Vol. 49. No.
	     10. ACM, 2014.
  [3] https://github.com/lemire/testingRNG/blob/master/source/splitmix64.h
  [4] Austin Appleby. Smhasher and murmurhash3 webpage, 2016.
      https://github.com/aappleby/smhasher/wiki
  [5] https://github.com/aappleby/smhasher/wiki/MurmurHash3
"""

# broccohash--
#	seed: a 64-bit seed for the hash function
#	v:    the 64-bit value to hash

def broccohash(seed,v):  # (both arguments are 64-bit values)
	seed ^= 0x3243F6A8885A308D     # (pi in base 16)
	u = v + seed
	u ^= u >> 30
	u = (u * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
	u ^= u >> 27
	u += seed >> 5
	u = (u * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
	u ^= u >> 31
	return u
