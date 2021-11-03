#!/usr/bin/env python3

from sys import stderr


# set_up_hash_function--
#	Import the appropriate hash function module and create the hasher

minimap2_hash = None
murmurhash3   = None
broccohash    = None

minimap2_hash_names = ["minimap2","minimap2hash","minimap2_hash"]
murmurhash3_names   = ["murmurhash3"]
broccohash_names    = ["brocco","broccohash","brocco_hash","Brocco","BroccoHash"]

def set_up_hash_function(hashType,hashSeed,kmerSize):
	global minimap2_hash,murmurhash3,broccohash
	if (hashType in minimap2_hash_names):
		if (minimap2_hash == None):
			try:
				from minimap2_hash import minimap2_hash
			except ImportError:
				print("WARNING: module minimap2_hash not found, using minimap2_hash_uncompiled",file=stderr)
				from minimap2_hash_uncompiled import minimap2_hash
		hashMask = (4**kmerSize)-1
		return lambda kmerBits : minimap2_hash(hashSeed,kmerBits,hashMask)
	elif (hashType in murmurhash3_names):
		if (murmurhash3 == None):
			try:
				from murmurhash3 import murmurhash3
			except ImportError:
				print("WARNING: module murmurhash3 not found, using murmurhash3_uncompiled",file=stderr)
				from murmurhash3_uncompiled import murmurhash3
		return lambda kmerBits : murmurhash3(hashSeed,kmerBits)
	elif (hashType in broccohash_names):
		if (broccohash == None):
			try:
				from broccohash import broccohash
			except ImportError:
				print("WARNING: module broccohash not found, using broccohash_uncompiled",file=stderr)
				from broccohash_uncompiled import broccohash
		return lambda kmerBits : broccohash(hashSeed,kmerBits)
	else:
		assert (False), "only minimap2 hash, murmurhash3, and broccohash are currently supported, not \"%s\"" % hashType
