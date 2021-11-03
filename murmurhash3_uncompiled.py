#!/usr/bin/env python3
"""
Wrapper for the 64-bit murmurhash3 function.

This is an uncompiled version to use if murmurhash3.pyx is unavailable or can't
be compiled.

The underlying hash function, murmurhash3(), was created by Austin Appleby.
The implementation here is derived from the mash source code (the function
MurmurHash3_x64_128 in common/murmur3.h). It is modified here as follows:
  - by hardwiring the length of the input to 64 bits (real murmurhash3 operates
    on an input buffer of any length)
  - by reduction of the output to 64 bits (instead of 128 bit)

References
  [1] MurmurHash, at wikipedia
      https://en.wikipedia.org/wiki/MurmurHash#MurmurHash3
  [2] MashMap, at github
      https://github.com/marbl/MashMap
  [3] mmh3, a Python wrapper for MurmurHash (MurmurHash3)
      https://github.com/hajimes/mmh3
      https://github.com/hajimes/mmh3/blob/master/mmh3module.cpp

Notes
  (1) Collisions *are* possible with this hash function.
  (2) The byte order is well-determined in this implementation, thus it will
      give the same results on all platforms. However, standard implementations
      of 128-bit MurmurHash3 will give different results on little-endian vs
      little-endian machines.
  (3) This implementation has been tested against version 2.5.1 of [3].
      Specifically, h1 and h2 are identical (by spot checking) to what [3]
      produces for
        mmh3.hash64(int.to_bytes(v,8,"little"),seed=s,signed=False)
      where v is a 64-bit value and s is a 32-bit seed. The implementation
      here accepts 64-bit seeds, while [3] appears to only use 32 bits of the
      seed.
"""

# murmurhash3--
#	seed: a 64-bit seed for the hash function
#	v:    the 64-bit value to hash

def murmurhash3(seed,v,which="xor"):  # (both arguments are 64-bit values)
	seed &= 0xFFFFFFFFFFFFFFFF
	h1 = seed
	h2 = seed

	c1 = 0x87C37B91114253D5 
	c2 = 0x4CF5AD432745937F

	# body
	# (normally murmurhash3 would process 128-bit blocks here but we have no
	# such blocks)

	# tail
	# (normally murmurhash3 would branch based on the number of bytes in
	# the final partial block but we always have exactly 8 bytes)

	#NA tail = pointer to bytes in v

	#NA k1 = 0
	#NA k2 = 0

	#NA k2 *= c2                                (not necessary because
	#NA k2  = ROTL64(k2,33)                      k2 will remain zero)
	#NA k2 *= c1
	#NA h2 ^= k2

	k1 = v & 0xFFFFFFFFFFFFFFFF							# k1 ^= tail[7] << 56
														# k1 ^= tail[6] << 48
														#  ...
														# k1 ^= tail[1] << 8
														# k1 ^= tail[0] << 0
	k1 = (k1 * c1) & 0xFFFFFFFFFFFFFFFF					# k1 *= c1
	k1 = ((k1 << 31) | (k1 >> 33)) & 0xFFFFFFFFFFFFFFFF	# k1 = ROTL64(k1,31)
	k1 = (k1 * c2) & 0xFFFFFFFFFFFFFFFF					# k1 *= c2
	h1 ^= k1

	# finalization

	# (murmurhash3 mixes in the length to prevent a certain type of collision;
	# in our case, since length is always the same, that type of collision
	# would never happen; so this is not really necessary but we do it anyway
	# so we can validate our output by comparison to other implementations)
	h1 ^= 8												# h1 ^= len
	h2 ^= 8												# h2 ^= len

	h1 = (h1 + h2) & 0xFFFFFFFFFFFFFFFF					# h1 += h2
	h2 = (h2 + h1) & 0xFFFFFFFFFFFFFFFF					# h2 += h1

	h1 = fmix64(h1)
	h2 = fmix64(h2)

	h1 = (h1 + h2) & 0xFFFFFFFFFFFFFFFF					# h1 += h2
	h2 = (h2 + h1) & 0xFFFFFFFFFFFFFFFF					# h2 += h1

	if (which == "h1"): return h1
	if (which == "h2"): return h2
	if (which == "+"):  return (h1+h2) & 0xFFFFFFFFFFFFFFFF
	# if (which == "xor"):                              # (we default to xor)
	return h1 ^ h2										# (real murmurhash3 returns 128 bits here instead)


def fmix64(k):  # (k is a 64-bit value)
	k ^= k >> 33
	k = (k * 0xFF51AFD7ED558CCD) & 0xFFFFFFFFFFFFFFFF	# k *= 0xFF51AFD7ED558CCD
	k ^= k >> 33
	k = (k * 0xC4CEB9FE1A85EC53) & 0xFFFFFFFFFFFFFFFF	# k *= 0xC4CEB9FE1A85EC53
	k ^= k >> 33
	return k

