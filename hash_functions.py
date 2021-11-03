#!/usr/bin/env python3

from sys import stderr


# set_up_hash_function--
#	Import the appropriate hash function module and create the hasher

minimap2_hash   = None
murmurhash3     = None
splitmix64_hash = None

minimap2_hash_names   = ["minimap2","minimap2hash","minimap2_hash"]
murmurhash3_names     = ["murmurhash3","MurmurHash3"]
splitmix64_hash_names = ["splitmix64","splitmix64hash","splitmix64_hash","SplitMix64","SplitMix64Hash"]

def set_up_hash_function(hashType,hashSeed,kmerSize):
	global minimap2_hash,murmurhash3,splitmix64_hash
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
	elif (hashType in splitmix64_hash_names):
		if (splitmix64_hash == None):
			try:
				from splitmix64_hash import splitmix64_hash
			except ImportError:
				print("WARNING: module splitmix64_hash not found, using splitmix64_hash_uncompiled",file=stderr)
				from splitmix64_hash_uncompiled import splitmix64_hash
		return lambda kmerBits : splitmix64_hash(hashSeed,kmerBits)
	else:
		assert (False), "only minimap2 hash, murmurhash3, and splitmix64_hash are currently supported, not \"%s\"" % hashType
