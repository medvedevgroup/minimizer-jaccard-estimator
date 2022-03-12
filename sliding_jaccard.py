#!/usr/bin/env python3

from sys                 import argv,stdin,stdout,stderr,exit
from math                import floor,ceil
from gzip                import open as gzip_open
from hash_functions      import set_up_hash_function, \
                                minimap2_hash_names, \
                                murmurhash3_names, \
                                splitmix64_hash_names
from winnowed_minimizers import winnowed_minimizers_linear

programName    = "sliding_jaccard"
programVersion = "0.5.1"


def usage(s=None):
	message = """
usage: cat <reference_fasta_file> | %s <query_fasta_file> [options]
  <query_fasta_file>       (required) file containing the sequence to compare
                           vs the <fasta_file>; this is expected to be 'short',
                           e.g. about 1Kbp long
  --window=<N>             (W=) minimizer window size (number of kmers in a
                           window)
                           (default is 100)
  --k=<N>                  (K=) kmer size
                           (default is 16)
  --canonical              consider reverse-complemented equivalent kmers to be
                           the same
                           (by default we consider such kmers as different)
  --hash=[<type>.]<seed>   type and seed for hash function; seed is an integer,
                           and "0x" prefix can be used to indicate hexadecimal;
                           type is either minimap2, murmurhash3, or splitmix64
                           (default is minimap2.0)
  --minimizers:local       minimizers are freshly recomputed for every
                           reference window     
                           (this is the default)
  --minimizers:global      minimizers are computed once for the reference, for
                           the entire sequence, and then intersected with each
                           window
  --distribution:<value>=<filename>
                           (cumulative) report the distribution of a particular
                           value; <value> is one of
                             J(Q,R)                   or true
                             J(Q,R;w)                 or winnowed
                             J(Q,R)-J(Q,R;w)          or true-winnowed
                             (J(Q,R)-J(Q,R;w))/J(Q,R) or (true-winnowed)/true
  --inhibit:details        don't report a table of J(Q,R) and J(Q,R;w) for each
                           position of the sliding window in the reference
                           (by default we report that, to stdout)
  --head=<number>          limit the number of input sequences
  --progress=<number>      periodically report how many sequences we've
                           processed
  --progress:scan=<number> periodically report how MUCH of the current
                           sequence we've processed
  --version                report this program's version number

Compute the Jaccard index between a query and every same-length of the
reference.""" \
% programName

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportDistribution
	global debug

	# parse the command line

	queryFilename        = None
	windowSize           = 100
	kmerSize             = 16
	canonical            = False
	hashType             = "minimap2"
	hashSeed             = 0
	globalMinimizers     = False
	reportDistribution   = None
	reportSlidingDetails = True
	sequenceLengthLimit  = None
	headLimit            = None
	reportProgress       = None
	reportScanProgress   = None
	debug                = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--help","--h","-help","-h"]):
			usage()
		elif (arg in ["--version","--v","--V","-version","-v","-V"]):
			exit("%s, version %s" % (programName,programVersion))
		elif (arg.startswith("--window=")) or (arg.startswith("W=")) or (arg.startswith("w=")):
			windowSize = int_with_unit(argVal)
			if (windowSize < 2):
				usage("window size has to be at least 2")
		elif (arg.startswith("--kmer=")) or (arg.startswith("K=")) or (arg.startswith("k=")):
			kmerSize = int(argVal)
			if (kmerSize < 2):
				usage("kmer size has to be at least 2")
			if (kmerSize > 32):
				usage("kmer size can't be more than 32 (because of hash function implementations)")
		elif (arg in ["--canonical","--canonicalize","--canon"]):
			canonical = True
		elif (arg.startswith("--hash=")) or (arg.startswith("--hasher=")):
			if (argVal in minimap2_hash_names):
				hashType = "minimap2"
				hashSeed = 0
			elif (argVal in murmurhash3_names):
				hashType = "murmurhash3"
				hashSeed = 0
			elif (argVal in splitmix64_hash_names):
				hashType = "splitmix64"
				hashSeed = 0
			else:
				if ("." in argVal):
					(hashType,hashSeed) = argVal.split(".",1)
				else:
					hashType = "minimap2"
					hashSeed = argVal
				if (hashSeed.startswith("0x")):
					hashSeed = int(hashSeed[2:],16)
				else:
					hashSeed = int(hashSeed)
		elif (arg == "--minimizers:local"):
			globalMinimizers = False
		elif (arg == "--minimizers:global"):
			globalMinimizers = True
		elif (arg.startswith("--distribution:J(Q,R)=")) \
		  or (arg.startswith("--distribution:true=")):
			if (reportDistribution == None): reportDistribution = {}
			reportDistribution["J(Q,R)"] = argVal.strip()
		elif (arg.startswith("--distribution:J(Q,R;w)=")) \
		  or (arg.startswith("--distribution:winnowed=")):
			if (reportDistribution == None): reportDistribution = {}
			reportDistribution["J(Q,R;w)"] = argVal.strip()
		elif (arg.startswith("--distribution:J(Q,R)-J(Q,R;w)=")) \
		  or (arg.startswith("--distribution:true-winnowed=")):
			if (reportDistribution == None): reportDistribution = {}
			reportDistribution["J(Q,R)-J(Q,R;w)"] = argVal.strip()
		elif (arg.startswith("--distribution:(J(Q,R)-J(Q,R;w))/J(Q,R)=")) \
		  or (arg.startswith("--distribution:(true-winnowed)/true=")):
			if (reportDistribution == None): reportDistribution = {}
			reportDistribution["(J(Q,R)-J(Q,R;w))/J(Q,R)"] = argVal.strip()
		elif (arg == "--inhibit:details"):
			reportSlidingDetails = False
		elif (arg.startswith("--maxlength=")):  # undocumented
			sequenceLengthLimit = int_with_unit(argVal)
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg.startswith("--progress:scan=")):
			reportScanProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (queryFilename == None):
			queryFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (queryFilename == None):
		usage("you must provide a query filename")

	reportSpecifics = ("specifics" in debug)
	if (reportSpecifics) and (globalMinimizers):
		exit ("--debug=specifics is only implemented for --minimizers:local")

	# set up the hash function and winnower

	hashFunc = set_up_hash_function(hashType,hashSeed,kmerSize)
	winnower = winnowed_minimizers_linear

	# read and preprocess the query

	if (queryFilename.endswith(".gz")) or (queryFilename.endswith(".gzip")):
		queryF = gzip_open(queryFilename,"rt")
	else:
		queryF = open(queryFilename,"rt")

	querySeq = None
	for (name,seq) in fasta_sequences(queryF):
		if (querySeq == None):
			(queryName,querySeq) = (name,seq.upper())
		else:
			assert (False), \
			       "\"%s\" contains more than one sequence" \
			     % queryFilename

	nonACGT = sum([1 for nuc in querySeq if (nuc not in "ACGT")])
	if (nonACGT > 0):
		print("WARNING: \"%s\" contains non-ACGT; kmerization will use a special hash value for these" % queryName,file=stderr)

	queryLen = len(querySeq)
	bufferLen = queryLen - (kmerSize-1)

	queryKmerProfile = kmer_profile(querySeq,kmerSize,canonical=canonical)

	if (reportScanProgress != None):
		print("hashing \"%s\"" % queryName,file=stderr)
	hashQ = hash_sequence(querySeq,kmerSize,hashFunc,canonical=canonical)
	if ("hashes" in debug):
		for (ix,h) in enumerate(hashQ):
			print("query h[%d] %016X" % (ix,h),file=stderr)

	if (reportScanProgress != None):
		print("winnowing \"%s\"" % queryName,file=stderr)
	queryMiniProfile = minimizer_profile(winnower(hashQ,windowSize))
	if ("minis" in debug):
		for (miniH,ix) in winnower(hashQ,windowSize):
			print("query mini %016X %d" % (miniH,ix),file=stderr)
		miniVals = [h for h in queryMiniProfile]
		miniVals.sort()
		for h in miniVals:
			if (queryMiniProfile[h] == 1):
				print("query profile %016X" % h,file=stderr)
			else:
				print("query profile %016Xx%d" % (h,queryMiniProfile[h]),file=stderr)

	if (reportSpecifics):
		queryMiniPositionToHash = {}
		for (miniH,ix) in winnower(hashQ,windowSize):
			queryMiniPositionToHash[ix] = miniH

	if ("bufferlength" in debug):
		print("queryLen=%d bufferLen=%d" % (queryLen,bufferLen),file=stderr)

	# process the reference sequences

	if (reportSlidingDetails):
		line = "#qName\tqLen" \
		     + "\trName\trStart\trEnd" \
		     + "\tw\tk" \
		     + "\tI(Q,R)\tU(Q,R)\tJ(Q,R)" \
		     + "\tI(Q,R;w)\tU(Q,R;w)\tJ(Q,R;w)"
		if (reportSpecifics): line += "\tminimizerLocations"
		print(line)

	seqNumber = 0
	for (name,seq) in fasta_sequences(stdin):
		seqNumber += 1
		if (headLimit != None) and (seqNumber > headLimit):
			print("limit of %d sequences reached" % headLimit,file=stderr)
			break
		if (reportProgress != None) \
		    and ((seqNumber == 1) or (seqNumber % reportProgress == 0)):
			print("processing sequence #%d: %s vs %s" % (seqNumber,name,nameB),file=stderr)

		seq = seq.upper()
		seqLen = len(seq)
		if (sequenceLengthLimit != None) and (seqLen >= sequenceLengthLimit):
			print("WARNING: truncating \"%s\" to %s bp (real length is %s bp)" \
			    % (name,"{:,}".format(sequenceLengthLimit),"{:,}".format(seqLen)),
			      file=stderr)
			seq = seq[:sequenceLengthLimit]
		if (seqLen < queryLen):
			print("WARNING: \"%s\" is shorter than the query" % name,file=stderr)
			continue

		nonACGT = sum([1 for nuc in seq if (nuc not in "ACGT")])
		if (nonACGT > 0):
			print("WARNING: \"%s\" contains non-ACGT; kmerization will use a special hash value for these" % name,file=stderr)

		if (reportDistribution != None):
			trueJaccardToCount        = {}    # J(Q,R)
			trueJaccardBucketSize     = 0.01
			winnowedJaccardToCount    = {}    # J(Q,R;w))
			winnowedJaccardBucketSize = 0.01
			jaccardErrorToCount       = {}    # J(Q,R)-J(Q,R;w))
			jaccardErrorBucketSize    = 0.01
			jaccardRelErrorToCount    = {}    # (J(Q,R)-J(Q,R;w))/J(Q,R)
			jaccardRelErrorBucketSize = 0.01

		# initialize the profiles

		kmerProfile = SlidingProfile(queryKmerProfile)
		miniProfile = SlidingProfile(queryMiniProfile)

		# initialize reference minimizer scanning

		if (reportScanProgress != None):
			print("hashing \"%s\"" % name,file=stderr)
		hashR = hash_sequence(seq,kmerSize,hashFunc,canonical=canonical)
		if ("hashes" in debug):
			for (ix,h) in enumerate(hashR):
				print("ref h[%d] %016X" % (ix,h),file=stderr)

		if (globalMinimizers):
			if (reportScanProgress != None):
				print("winnowing \"%s\"" % name,file=stderr)
			refMinis = list(winnower(hashR,windowSize))
			if ("minis" in debug):
				for (miniH,seqIx) in refMinis:
					print("ref mini %016X %d" % (miniH,seqIx),file=stderr)

			assert (refMinis != [])
			oldMiniIx   = 0
			nextMiniIx  = 0
			oldMiniPos  = refMinis[oldMiniIx][1] if (oldMiniIx < len(refMinis)) \
						  else None
			nextMiniPos = refMinis[nextMiniIx][1] if (nextMiniIx < len(refMinis)) \
						  else None

			if ("minis" in debug):
				if (nextMiniPos != None):
					(h,seqIx) = refMinis[nextMiniIx]
					print("ref mini next:   [%d] hashR[%d]=%016X" % (nextMiniIx,seqIx,h),file=stderr)
				else:
					print("ref mini next:   (none)",file=stderr)
				if (oldMiniPos != None):
					(h,seqIx) = refMinis[oldMiniIx]
					print("ref mini old:    [%d] hashR[%d]=%016X" % (oldMiniIx,seqIx,h),file=stderr)
				else:
					print("ref mini old:    (none)",file=stderr)

		# process the reference, kmer-by-kmer

		kmerBuffer = []
		kmerPos = -1
		for kmer in kmerize(seq,kmerSize,canonical=canonical):
			kmerBuffer += [kmer]
			kmerPos += 1
			refStart = kmerPos+1 - bufferLen
			refEnd   = kmerPos + kmerSize

			if (reportScanProgress != None) \
			    and ((kmerPos == 1) or (kmerPos % reportScanProgress == 0)):
				print("scanning position %s" % ("{:,}".format(kmerPos)),file=stderr)

			if ("minis" in debug):
				print("kmerPos=%d hashR[%d]=%016X" \
				    % (kmerPos,kmerPos,hashR[kmerPos]),
				      file=stderr)

			# (for global minimizers) if we have reached a new minimizer; incorporate it

			if (globalMinimizers):
				if (nextMiniPos != None) and (kmerPos == nextMiniPos):
					miniH = refMinis[nextMiniIx][0]
					if ("minis" in debug):
						print("ref mini add:    [%d] hashR[%d]=%016X" % (nextMiniIx,nextMiniPos,miniH),file=stderr)

					miniProfile.add(miniH)

					nextMiniIx += 1
					nextMiniPos = refMinis[nextMiniIx][1] if (nextMiniIx < len(refMinis)) \
								  else None

					if ("minis" in debug):
						print("mini.nI=%d mini.nU=%d" % (miniProfile.nI,miniProfile.nU),file=stderr)
						if (nextMiniPos != None):
							(h,seqIx) = refMinis[nextMiniIx]
							print("ref mini next:   [%d] hashR[%d]=%016X" % (nextMiniIx,seqIx,h),file=stderr)
						else:
							print("ref mini next:   (none)",file=stderr)

				# if a minimizer is now out of the window; dis-incorporate it

				if (oldMiniPos != None) and (kmerPos-bufferLen == oldMiniPos):
					miniH = refMinis[oldMiniIx][0]
					if ("minis" in debug):
						print("ref mini remove: [%d] hashR[%d]=%016X" % (oldMiniIx,oldMiniPos,miniH),file=stderr)

					miniProfile.remove(miniH)

					oldMiniIx += 1
					oldMiniPos = refMinis[oldMiniIx][1] if (oldMiniIx < len(refMinis)) \
								 else None

					if ("minis" in debug):
						print("mini.nI=%d mini.nU=%d" % (miniProfile.nI,miniProfile.nU),file=stderr)
						if (oldMiniPos != None):
							(h,seqIx) = refMinis[oldMiniIx]
							print("ref mini old:    [%d] hashR[%d]=%016X" % (oldMiniIx,seqIx,h),file=stderr)
						else:
							print("ref mini old:    (none)",file=stderr)

			# incorporate the latest kmer into the profile

			kmerProfile.add(kmer)

			# if we don't have a full kmer buffer to compute from, go back and
			# get the next kmer

			if (len(kmerBuffer) < bufferLen):
				continue   

			# (for local minimizers) compute minimizers for this window and
			# incorporate them (*just* them) into the profile

			if (not globalMinimizers):
				windowStart = kmerPos+1-bufferLen
				if (reportScanProgress != None) \
				    and ((kmerPos == 1) or (kmerPos % reportScanProgress == 0)):
					print("winnowing \"%s\"[%d:%d]" \
					    % (name,windowStart,windowStart+bufferLen),
					      file=stderr)
				miniProfile.reset()
				miniPositionToHash = {}   # only used if reportSpecifics
				for (miniH,seqIx) in winnower(hashR[windowStart:windowStart+bufferLen],windowSize):
					if ("minis" in debug):
						print("ref[%d:%d] @%d mini %016X %d" \
						    % (windowStart,windowStart+bufferLen,seqIx,miniH,windowStart+seqIx),
						       file=stderr)
					miniProfile.add(miniH)
					miniPositionToHash[seqIx] = miniH

			# if the buffer is over-full, dis-incorporate the oldest kmer
			
			if (len(kmerBuffer) > bufferLen):
				kmerProfile.remove(kmerBuffer.pop(0))

			# report the jaccard of the profile -OR- collect the distribution

			kmerJaccard = kmerProfile.jaccard()
			miniJaccard = miniProfile.jaccard()

			specifics = None
			if (reportSpecifics):
				# build a string that shows the positions of all minimizers
				#  Q => minimizer in query (and also in reference sliding window)
				#  q => minimizer in query (but not in reference sliding window)
				#  S => minimizer in reference sliding window (and also in query)
				#  s => minimizer in reference sliding window (but not in query)
				#  * => minimizer at this position in both query and reference sliding window
				#  | => separator dividing off windowLength-1 positions at each end
				specifics = ["-"] * bufferLen
				for seqIx in queryMiniPositionToHash:
					h = queryMiniPositionToHash[seqIx]
					if (h in miniProfile):
						specifics[seqIx] = "Q"
					else:
						specifics[seqIx] = "q"
				for seqIx in miniPositionToHash:
					h = miniPositionToHash[seqIx]
					if (specifics[seqIx] != "-"):
						specifics[seqIx] = "*"
					elif (h in queryMiniProfile):
						specifics[seqIx] = "S"
					else:
						specifics[seqIx] = "s"
				specifics = "".join(specifics[:(windowSize-1)]) \
				          + "|" + "".join(specifics[(windowSize-1):-(windowSize-1)]) \
				          + "|" + "".join(specifics[-(windowSize-1):])

			if (reportSlidingDetails):
				line = ("%s\t%d"
				      + "\t%s\t%d\t%d"
				      + "\t%d\t%d"
				      + "\t%d\t%d\t%.6f"
				      + "\t%d\t%d\t%.6f") \
				     % (queryName,queryLen,
				        name,refStart,refEnd,
				        windowSize,kmerSize,
				        kmerProfile.nI,kmerProfile.nU,kmerJaccard,
				        miniProfile.nI,miniProfile.nU,miniJaccard)
				if (specifics != None): line += "\t" + specifics
				print(line)

			if (reportDistribution != None):
				bucket = value_to_bucket(kmerJaccard,trueJaccardBucketSize)
				if (bucket not in trueJaccardToCount): trueJaccardToCount[bucket] =  1
				else:                                  trueJaccardToCount[bucket] += 1

				bucket = value_to_bucket(miniJaccard,winnowedJaccardBucketSize)
				if (bucket not in winnowedJaccardToCount): winnowedJaccardToCount[bucket] =  1
				else:                                      winnowedJaccardToCount[bucket] += 1

				jaccardError = kmerJaccard - miniJaccard
				bucket = value_to_bucket(jaccardError,jaccardErrorBucketSize)
				if (bucket not in jaccardErrorToCount): jaccardErrorToCount[bucket] =  1
				else:                                   jaccardErrorToCount[bucket] += 1

				# $$$ is this the right way to handle a zero denominator?
				if (kmerJaccard != 0):
					jaccardRelError = jaccardError / kmerJaccard
					bucket = value_to_bucket(jaccardRelError,jaccardRelErrorBucketSize)
					if (bucket not in jaccardRelErrorToCount): jaccardRelErrorToCount[bucket] =  1
					else:                                      jaccardRelErrorToCount[bucket] += 1

		if (reportDistribution != None):
			report_distribution("J(Q,R)",                  trueJaccardToCount,    trueJaccardBucketSize,    queryName,queryLen,name,seqLen)
			report_distribution("J(Q,R;w)",                winnowedJaccardToCount,winnowedJaccardBucketSize,queryName,queryLen,name,seqLen)
			report_distribution("J(Q,R)-J(Q,R;w)",         jaccardErrorToCount,   jaccardErrorBucketSize,   queryName,queryLen,name,seqLen)
			report_distribution("(J(Q,R)-J(Q,R;w))/J(Q,R)",jaccardRelErrorToCount,jaccardRelErrorBucketSize,queryName,queryLen,name,seqLen)

	if (reportDistribution != None):
		for variable in reportDistribution:
			f = reportDistribution[variable]
			if (type(f) != str):
				f.close()


# SlidingProfile--
#	Class to track the intersection/union between a static profile and a
#	profile for a sliding window. Items in the profile can be any basic type,
#	for example kmers or hash values.

class SlidingProfile(object):
	# profile is an iteratable that represents a set of elements of some type;
	# for example, it could be a set of kmers or a list of minimizer hash
	# values

	def __init__(self,profile):
		self.staticProfile  = set(profile)
		self.reset()

	def reset(self):
		self.rollingProfile = {}                          # maps item to count
		self.nI = 0                                       # intersection
		self.nU = sum([1 for item in self.staticProfile]) # union

	def add(self,item):
		# incorporate an item into the rolling profile
		if (item not in self.rollingProfile): self.rollingProfile[item] =  1
		else:                                 self.rollingProfile[item] += 1

		if (item not in self.staticProfile):
			if (self.rollingProfile[item] == 1): self.nU += 1
		else: # if (item in self.staticProfile):
			if (self.rollingProfile[item] == 1): self.nI += 1

	def remove(self,item):
		# dis-incorporate an item from the rolling profile
		self.rollingProfile[item] -= 1
		if (item not in self.staticProfile):
			if (self.rollingProfile[item] == 0): self.nU -= 1
		else: # if (item in self.staticProfile):
			if (self.rollingProfile[item] == 0): self.nI -= 1

	def __contains__(self,item):
		return (item in self.rollingProfile) and (self.rollingProfile[item] > 0)

	def __getitem__(self,item):
		if (item in self.rollingProfile): return self.rollingProfile[item]
		else:                             return 0

	def jaccard(self):
		return 0.0 if (self.nU == None) \
		       else self.nI/self.nU                       # (non-truncating divide)


# kmer_profile--
#	Return a dict, mapping a kmer to number of occurrences

def kmer_profile(seq,kmerSize,canonical=False):
	seqLen = len(seq)
	if (canonical):
		revSeq = reverse_complement(seq)

	profile = {}
	for ix in range(kmerSize,seqLen+1):
		kmer = seq[ix-kmerSize:ix]
		if (canonical):
			revKmer = revSeq[seqLen-ix:seqLen-ix+kmerSize]
			kmer = min(kmer,revKmer)

		if (kmer not in profile): profile[kmer] =  1
		else:                     profile[kmer] += 1

	return profile

# kmerize--
#	Yield the kmers in a sequence, one-by-one

def kmerize(seq,kmerSize,canonical=False):
	seqLen = len(seq)
	if (canonical):
		revSeq = reverse_complement(seq)

	for ix in range(kmerSize,seqLen+1):
		kmer = seq[ix-kmerSize:ix]
		if (canonical):
			revKmer = revSeq[seqLen-ix:seqLen-ix+kmerSize]
			kmer = min(kmer,revKmer)

		yield kmer


# hash_sequence--
#	Convert a sequence to hashed kmers

ntToBits = {"A":0 , "C":1, "G":2, "T":3}
hashOfBadKmer = 0xFFFFFFFFFFFFFFFF   # assumes hash range is 64 bits

def hash_sequence(seq,kmerSize,hashFunc,canonical=False):
	seqLen = len(seq)

	if (canonical):
		revSeq = reverse_complement(seq)

	hashSeq = []
	for ix in range(kmerSize,seqLen+1):
		kmer = seq[ix-kmerSize:ix]
		if (canonical):
			revKmer = revSeq[seqLen-ix:seqLen-ix+kmerSize]
			kmer = min(kmer,revKmer)

		try:
			kmerBits = 0
			for nt in kmer:
				kmerBits = (kmerBits<<2) + ntToBits[nt]
			hashSeq += [hashFunc(kmerBits)]
		except KeyError:
			# the kmer contains a non-ACGT
			hashSeq += [hashOfBadKmer]

	return tuple(hashSeq)


# minimizer_profile--
#	Return a dict, mapping a hash value to number of occurrences

def minimizer_profile(minimizers):
	profile = {}
	for (h,ix) in minimizers:
		if (h not in profile): profile[h] =  1
		else:                  profile[h] += 1

	return profile


# fasta_sequences--

def fasta_sequences(f):
	seqName = None
	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seq))
			(seqName,seq) = (line[1:].strip(),[])
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seq += [line]

	if (seqName != None):
		yield (seqName,"".join(seq))


# report_distribution--

def report_distribution(variable,variableToCount,bucketSize,queryName,queryLen,refName,refLen):
	global reportDistribution

	if (variable in reportDistribution):
		f = reportDistribution[variable]
		if (type(f) == str):
			filename = f
			if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
				f = gzip_open(filename,"wt")
			else:
				f = open(filename,"wt")
			reportDistribution[variable] = f

			print("#qName\tqLen\trName\trLen\t%s\tcount(%s)"
			    % (variable,variable),
				  file=f)

		buckets = [bucket for bucket in variableToCount]
		buckets.sort()
		for bucket in buckets:
			bucketMiddle = bucket_to_value(bucket,bucketSize)
			print(("%s\t%d\t%s\t%d\t%.6f\t%d"
				% (queryName,queryLen,refName,refLen,
				   bucketMiddle,variableToCount[bucket])),
				  file=f)

		f.close()


# value_to_bucket, bucket_to_value--
#	For binning distributions, convert a value to the bucket in which it
#	should be counted, or convert a bucket to the value at the bucket's center

def value_to_bucket(v,bucketSize):
	return int(floor((v + (bucketSize/2)) / bucketSize))

def bucket_to_value(bucket,bucketSize):
	return (bucketSize*bucket)


# reverse_complement--

complementMap = str.maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                              "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
