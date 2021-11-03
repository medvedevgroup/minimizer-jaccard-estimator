#!/usr/bin/env python3
"""
Compute the Jaccard correction for pairs of sequences, relating to winnowed
minimizers.
"""

from sys                 import argv,stdin,stdout,stderr,exit
from math                import ceil
from gzip                import open as gzip_open
from random              import seed as random_seed,randint
from itertools           import product as cartesian_product
from collections         import defaultdict
from time                import strftime
from hash_functions      import set_up_hash_function, \
                                minimap2_hash_names, \
                                murmurhash3_names, \
                                broccohash_names
from winnowed_minimizers import winnowed_minimizers_linear

programName    = "jaccard_correction_test"
programVersion = "0.2.3"


def usage(s=None):
	message = """
usage: cat <fasta_file_1> | %s [<fasta_file_2>] [options]
  <fasta_file_2>          (optional) file to compare to (see notes below)
  --window=<N>            (W=) minimizer window size (number of kmers in a
                          window)
                          (default is 100)
  --k=<N>                 (K=) kmer size
                          (default is 16)
  --canonical             consider reverse-complemented equivalent kmers to be
                          the same
                          (by default we consider such kmers as different)
  --replicates=<N>        perform N replicates for each sequence pair, each
                          using a different hash function; see note below on
                          how hash functions are seeded
                          (default is 1 replicate)
  --hash=[<type>.]<seed>  type and seed for hash function; seed is an integer,
                          and "0x" prefix can be used to indicate hexadecimal;
                          type is either minimap2, murmurhash3, or broccohash
                          (default is minimap2.0)
  --prng=<string>         set seed for PRNG; NOT to be confused with the hash
                          seed
  --inhibit:correction    don't compute any of the statistics related to the
                          jaccard correction
  --head=<number>         limit the number of input sequence pairs
  --progress=<number>     periodically report how many sequence pairs we've
                          processed
  --version               report this program's version number

Compute the Jaccard correction for pairs of sequences, relating to winnowed
minimizers.

If fasta_file_2 is provided, the first sequence in fasta_file_1 is compared to
the first from fasta_file_2, the second is compare to the second, and so on.

If fasta_file_2 is not provided, the first sequence in fasta_file_1 is compared
to the second sequence, the third is compared to the fourth, and so on.

If there is only one replicate, it is seeded with the user's --hash=<seed>
setting. When there is more than one replicate, a PRNG is used to select N
distinct seeds. The first of these is the same as the user's --hash=<seed>
setting. The replicates for each sequence pair use the same sequence of hash
seeds.""" \
% programName

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	fasta2Name        = None
	windowSize        = 100
	kmerSize          = 16
	canonical         = False
	numReplicates     = 1
	hashType          = "minimap2"
	hashSeed          = 0
	prngSeed          = None
	computeCorrection = True
	headLimit         = None
	reportProgress    = None
	debug             = []

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
		elif (arg in ["--noncanonical","--nocanonicalize","--noncanon","--nocanon"]):
			canonical = False
		elif (arg.startswith("--replicates=")) or (arg.startswith("--reps=")):
			numReplicates = int(argVal)
			if (numReplicates < 1):
				usage("number of replicates has to be at strictly positive")
		elif (arg.startswith("--hash=")) or (arg.startswith("--hasher=")):
			if (argVal in minimap2_hash_names):
				hashType = "minimap2"
				hashSeed = 0
			elif (argVal in murmurhash3_names):
				hashType = "murmurhash3"
				hashSeed = 0
			elif (argVal in broccohash_names):
				hashType = "broccohash"
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
		elif (arg.startswith("--prng=")):
			if (argVal.lower() == "none"):
				prngSeed = None
			else:
				prngSeed = argVal
				prngSeed = prngSeed.replace("{date}",strftime("%d/%m/%Y"))
				prngSeed = prngSeed.replace("{time}",strftime("%I/%M/%S"))
		elif (arg in ["--inhibit:correction","--nocorrection"]):
			computeCorrection = False
		elif (arg in ["--noinhibit:correction","--correction"]):
			computeCorrection = True
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (fasta2Name == None):
			fasta2Name = arg
		else:
			usage("unrecognized option: %s" % arg)

	# set up the hash function(s) and winnower

	if (prngSeed != None):
		random_seed(prngSeed)

	if (numReplicates == 1):
		hashSeeds = [hashSeed]
		hashFuncs = [set_up_hash_function(hashType,hashSeed,kmerSize)]
	elif (numReplicates > 1):
		hashSeeds = []
		hashXor = None
		for _ in range(numReplicates):
			hashSeeds += [randint(0,(1<<64)-1)]
			if (hashXor == None):
				hashXor = hashSeeds[0] ^ hashSeed
				hashSeeds[0] = hashSeed
			else:
				hashSeeds[-1] ^= hashSeed
		hashFuncs = [set_up_hash_function(hashType,hashSeed,kmerSize) 
		             for hashSeed in hashSeeds]

	winnower = winnowed_minimizers_linear

	# open the second fasta file, if we have one

	fasta2F = None
	if (fasta2Name != None):
		if (fasta2Name.endswith(".gz")) or (fasta2Name.endswith(".gzip")):
			fasta2F = gzip_open(fasta2Name,"rt")
		else:
			fasta2F = open(fasta2Name,"rt")

	# process the fasta sequence pairs
	# $$$ "tau" was used in early versions of the manuscript, but it has been
	#     since changed to D
	# $$$ "alpha" as used in early versions of the manuscript, but it has been
	#     since changed to g

	showAllReplicates = ("replicates" in debug)

	header =  ["name1","name2"]
	if (numReplicates > 1):
		if (showAllReplicates): header += ["rep","hash.seed"]
		else:                   header += ["replicates"]
	header += ["W","K","length.nt","|a|"]
	header += ["I(A,B)","U(A,B)","J(A,B)"]
	if ("density" in debug): header += ["2L/(w+1)","|A.min(w)|","|B.min(w)|"]
	header += ["I(A,B;w)","U(A,B;w)","J(A,B;w)"]
	header += ["tau(A,B;w)","J(A,B)-eps","CI(A,B;w)","Bias(A,B;w)"]
	header += ["J(A,B;w)-J(A,B)","I(A,B;w)-CI(A,B;w)"]

	if ("checkduplicates" not in debug):
		print("#%s" % "\t".join(header))

	pairNumber = 0
	for (nameA,seqA,nameB,seqB) in fasta_pairs(stdin,fasta2F,fasta2Name):
		pairNumber += 1
		if (headLimit != None) and (pairNumber > headLimit):
			print("limit of %d pairs reached" % headLimit,file=stderr)
			break
		if (reportProgress != None) \
		    and ((pairNumber == 1) or (pairNumber % reportProgress == 0)):
			if (reportProgress == 1) and (numReplicates > 1):
				print("processing pair #%d: %s vs %s (replicate 1)" % (pairNumber,nameA,nameB),file=stderr)
			else:
				print("processing pair #%d: %s vs %s" % (pairNumber,nameA,nameB),file=stderr)

		assert (len(seqA) == len(seqB)), \
		       "%s and %s have different lengths (%d and %d)" \
		     % (nameA,nameB,len(seqA),len(seqB))

		kd = jaccard_by_kmers(seqA,seqB,kmerSize,canonical=canonical)
		# Sp^k(A) intersect Sp^k(B) is kd.nIntersection
		# Sp^k(A) union Sp^k(B) is kd.nUnion
		# J(A,B) is kd.jaccard

		nIntersectionSum = nUnionSum = 0
		jaccardSum = 0.0

		for replicateNum in range(numReplicates):
			if (reportProgress == 1) and (numReplicates > 1) and (replicateNum > 0):
				print("processing pair #%d: %s vs %s (replicate %d)"
				    % (pairNumber,nameA,nameB,1+replicateNum),file=stderr)

			hashFunc = hashFuncs[replicateNum]
			hashSeed = hashSeeds[replicateNum]

			hashA = hash_sequence(seqA,kmerSize,hashFunc,canonical=canonical)
			hashB = hash_sequence(seqB,kmerSize,hashFunc,canonical=canonical)
			md = jaccard_by_minimizers(hashA,hashB,windowSize,winnower)
			# I(A,B;w) is md.nIntersection
			# U(A,B;w) is md.nUnion
			# J(A,B;w) is md.jaccard

			nIntersectionSum += md.nIntersection
			nUnionSum        += md.nUnion
			jaccardSum       += md.jaccard

			if ("kmerseqs" in debug):
				dbg_report_kmer_seqs(seqA,seqB,kmerSize,canonical=canonical)
			if ("hashes" in debug):
				dbg_report_hashes(hashA,hashB,kmerSize)

			if ("checkduplicates" in debug):
				try:
					cd = jaccard_correction(nameA,hashA,nameB,hashB,windowSize)
					print("\"%s\" and \"%s\" are duplicate free" % (nameA,nameB),
					      file=stderr)
					# $$$ we also want to test that all matches are in the same
					#     .. position in both sequences
				except ValueError:
					pass
				continue

			if (computeCorrection) and (replicateNum == 0):
				# assuming that the hash function is invertible, this is only
				# needed for the first replicate (since it will be the same
				# result for all replicates)
				cd = jaccard_correction(nameA,hashA,nameB,hashB,windowSize)
				# tau is cd.tau
				# J(A,B)-epsilon is cd.jaccardByTau
				# C_I(A,B;w) is cd.CI
				# \mathcal B Bias(A,B;w) is cd.bias

				if ("CI" in debug):
					dbg_report_CI(cd.dbgCI)

			if (numReplicates > 1) and (not showAllReplicates):
				if (replicateNum < numReplicates-1):
					# (line printed for final replicate only)
					continue
				else:
					# (for final replicate we'll print these averages)
					nIntersectionAvg = nIntersectionSum / numReplicates
					nUnionAvg        = nUnionSum        / numReplicates
					jaccardAvg       = jaccardSum       / numReplicates

			line =  ["%s\t%s"               % (nameA,nameB)]
			if (numReplicates > 1):
				if (showAllReplicates): line += ["%d\t0x%016X" % (1+replicateNum,hashSeed)]
				else:                   line += ["%d"          % numReplicates]
			line += ["%d\t%d\t%d\t%d"       % (windowSize,kmerSize,len(seqA),len(hashA))]
			line += ["%d\t%d\t%.6f"         % (kd.nIntersection,kd.nUnion,kd.jaccard)]
			if ("density" in debug):
				line += ["%.3f\t%d\t%d"     % (2*len(hashA)/(windowSize+1),md.nMinimizersA,md.nMinimizersB)]
			if (numReplicates == 1) or (showAllReplicates):
				line += ["%d\t%d\t%.6f"     % (md.nIntersection,md.nUnion,md.jaccard)]
			elif (replicateNum == numReplicates-1):
				line += ["%.3f\t%.3f\t%.6f" % (nIntersectionAvg,nUnionAvg,jaccardAvg)]
			if (computeCorrection):
				line += ["%d\t%.6f\t%.6f\t%.6f" % (cd.tau,cd.jaccardByTau,cd.CI,cd.bias)]
			else:
				line += ["NA\tNA\tNA\tNA"]
			if (numReplicates == 1) or (showAllReplicates):
				if (computeCorrection):
					line += ["%.6f\t%.6f"   % (md.jaccard-kd.jaccard,md.nIntersection-cd.CI)]
				else:
					line += ["%.6f\tNA"     % (md.jaccard-kd.jaccard)]
			elif (replicateNum == numReplicates-1):
				if (computeCorrection):
					line += ["%.6f\t%.6f"   % (jaccardAvg-kd.jaccard,nIntersectionAvg-cd.CI)]
				else:
					line += ["%.6f\tNA"     % (jaccardAvg-kd.jaccard)]
			print("\t".join(line))

	if (fasta2F != None):
		fasta2F.close()


# jaccard_correction--
#	Compute the correction for jaccard, according to our manuscript

class CorrectionData: pass

def jaccard_correction(nameA,hashA,nameB,hashB,windowSize):

	# make sure the hash sequences contain no duplicates; failure would be do
	# either to having duplicate kmers in the nucleotide sequence, or using a
	# non-invertible hash such as murmurhash3

	try:
		duplicates = []
		aPosToBPos = shared_hash_positions(hashA,hashB,duplicates=duplicates)
	except ValueError:
		(which,h,pos1,pos2) = duplicates[0]
		if (which == "A"):
			print("\"%s\" contains duplicate hashed kmers; positions %d and %d both hash to 0x%016X" \
			    % (nameA,pos1,pos2,h),
			      file=stderr)
		else: # if (which == "B"):
			print("\"%s\" contains duplicate hashed kmers; positions %d and %d both hash to 0x%016X" \
			    % (nameB,pos1,pos2,h),
			      file=stderr)
		raise

	try:
		duplicates = []
		bPosToAPos = shared_hash_positions(hashB,hashA,duplicates=duplicates)
	except ValueError:
		(which,h,pos1,pos2) = duplicates[0]
		if (which == "A"):
			print("\"%s\" contains duplicate hashed kmers; positions %d and %d both hash to 0x%016X" \
			    % (nameA,pos1,pos2),
			      file=stderr)
		else: # if (which == "B"):
			print("\"%s\" contains duplicate hashed kmers; positions %d and %d both hash to 0x%016X" \
			    % (nameB,pos1,pos2),
			      file=stderr)
		raise

	# compute configurations

	shared = window_shared_counts(hashA,hashB,windowSize,aPosToBPos)
	# S(i,j,w) is shared[(i,j)], but only when S(i,j,w)>0

	if ("shared" in debug):
		dbg_report_shared(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize)

	nConfigurations = configuration_counts(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize)

	L = len(hashA)
	w = windowSize

	if ("CI" in debug):
		dbgCI = []

	tau = 0
	for s in range(windowSize+1):
		if (s not in nConfigurations): continue
		for (cal,cbl) in cartesian_product([0,1,2],[0,1,2]):
			tau += nConfigurations[s][(cal,0,cbl,0)]
	jTau = float(tau) / (2*L - tau)                  # J(A,B)-epsilon is jTau

	CI = 0                                           # C_I(A,B;w) is CI
	for s in range(windowSize+1):                    # (min is 0, max is w)
		f0 = 1.0 / (2*w-s)                           # f_0(s) is f0
		f1 = 1.0 / ((2*w-s)*(2*w-s+1))               # f_1(s) is f1
		f2 = 2.0 / ((2*w-s)*(2*w-s+1)*(2*w-s+2))     # f_2(s) is f2
		alpha = {}                                   # alpha_t(s) is alpha[t]
		for t in [0,1,2]:                            # (note that alpha_1(s) is not needed)
			alpha[t] = 0
			for (car,cbr) in cartesian_product([0,1,2],[0,1,2]):
				alpha[t] += (s-1)*nConfigurations[s][(t,car,t,cbr)]
			alpha[t] -= nConfigurations[s][(t,1,t,1)]
			alpha[t] += nConfigurations[s][(t,2,t,2)]
		if ("alpha" in debug):
			dbg_report_alpha(alpha,s)

		CI += f1*alpha[0] + f2*alpha[2]
		for (cal,cbl) in cartesian_product([0,1,2],[0,1,2]):
			CI += f0*nConfigurations[s][(cal,0,cbl,0)]
		for (cal,car,cbl) in cartesian_product([0,2],[0,1,2],[0,1,2]):
			CI += f1*nConfigurations[s][(cal,car,cbl,1)]
		for (cal,cbl,cbr) in cartesian_product([0,1,2],[0,2],[0,1,2]):
			CI += f1*nConfigurations[s][(cal,1,cbl,cbr)]

		if ("CI" in debug):
			deltaCI = CI if (s==0) else (CI-prevCI)
			prevCI = CI
			parts = []
			parts += [(f1,alpha[0])]
			parts += [(f2,alpha[2])]
			parts += [(f0,sum([nConfigurations[s][(cal,0,cbl,0)] for (cal,cbl) in cartesian_product([0,1,2],[0,1,2])]))]
			parts += [(f1,sum([nConfigurations[s][(cal,car,cbl,1)] for (cal,car,cbl) in cartesian_product([0,2],[0,1,2],[0,1,2])]))]
			parts += [(f1,sum([nConfigurations[s][(cal,1,cbl,cbr)] for (cal,cbl,cbr) in cartesian_product([0,1,2],[0,2],[0,1,2])]))]
			dbgCI += [(s,deltaCI,parts)]

	bias = (float((w+1)*CI)/(4*L-(w+1)*CI)) - (float(tau) / (2*L - tau))
	# \mathcal B Bias(A,B;w) is bias

	cd = CorrectionData()
	cd.tau          = tau
	cd.jaccardByTau = jTau
	cd.CI           = CI
	cd.bias         = bias
	if ("CI" in debug): cd.dbgCI = dbgCI
	return cd


# window_shared_counts
#	return S s.t. S[(i,j)] is the count of shared hash values in A[i..i+w-1]
#	and B[j..j+w-1]. S is only populated for (i,j) that lead to non-zero counts.

def window_shared_counts(hashA,hashB,windowSize,aPosToBPos):
	w = windowSize
	shared = {}
	for aPos in aPosToBPos:
		bPos = aPosToBPos[aPos]
		for i in range(aPos-(w-1),aPos+1):   # min is aPos-(w-1), max is aPos
			if (i < 0): continue
			if (i+w > len(hashA)): continue
			for j in range(bPos-(w-1),bPos+1):
				if (j < 0): continue
				if (j+w > len(hashB)): continue
				if ((i,j) not in shared): shared[(i,j)] =  1
				else:                     shared[(i,j)] += 1

	return shared


# configuration_counts--
#	Count the occurences of 'configurations' between the two sequences. We have
#	two implementations -- all_configurations() which is easy to prove correct,
#	and configurations() which is much faster.
#
# The vast majority of configs are [2,2,2,2]0. configurations() avoids
# reporting most of those, but computes all the others. The count for the
# unreported [2,2,2,2]0 is then deduced from the difference in the number of
# configurations reported and the expected number, configurations_count().
# Note that configurations() may report *some* [2,2,2,2]0.
#
# Hypothesis:
#	All non-[2,2,2,2]0 configurations fall into this set:
#	union over i',j' s.t. a[i']=b[j'] of {i,j s.t. i<=i'-w<=i and j<=j'-w<=j}
#
# The config at (i,j) is *not* [2,2,2,2]0 if any of these are true
#	s!=0     S(i+1,j+1,w) != 0           (i+1,j+1) in shared
#	cal!=2   a[i] is in b[j..j+w]        j <= aPosToBPos[i] <= j+w
#	car!=2   a[i+w] is in b[j+1..j+w]    j+1 <= aPosToBPos[i+w] <= j+w
#	cbl!=2   b[j] is in a[i..i+w]        i <= bPosToAPos[j] <= i+w
#	cbr!=2   b[j+w] is in a[i+1..i+w]    i+1 <= bPosToAPos[j+w] <= i+w

def configuration_counts(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize):
	configurations_reporter = all_configurations if ("allconfigs" in debug) \
	                     else configurations

	nExpectedConfigurations = configurations_count(hashA,hashB,windowSize)

	nConfigurations = {}
	for s in range(windowSize+1):        # min is 0, max is w
		nConfigurations[s] = defaultdict(int)

	ijToConfig = {}
	for (s,cal,car,cbl,cbr,i,j) in configurations_reporter(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize):
		nConfigurations[s][(cal,car,cbl,cbr)] += 1
		if ("diagconfigs" in debug):
			ijToConfig[(i,j)] = (s,cal,car,cbl,cbr)
		elif ("configurations" in debug):
			dbg_report_configuration(i,j,s,cal,car,cbl,cbr)

	if ("diagconfigs" in debug):
		w = windowSize
		maxDelta = len(hashA)-(w+1)
		for delta in range(-maxDelta,maxDelta+1):
			for i in range(len(hashA)-w):        # min is 0, max is |a|-(w+1)
				j = i + delta
				if (i,j) not in ijToConfig: continue
				(s,cal,car,cbl,cbr) = ijToConfig[(i,j)]
				dbg_report_configuration(i,j,s,cal,car,cbl,cbr)

	configurationsReported = 0
	for s in nConfigurations:
		for (cal,car,cbl,cbr) in nConfigurations[s]:
			configurationsReported += nConfigurations[s][(cal,car,cbl,cbr)]
	if ("allconfigs" in debug):
		assert (configurationsReported == nExpectedConfigurations)
	else:
		assert (configurationsReported <= nExpectedConfigurations)
		(s,cal,car,cbl,cbr) = unreportedConfiguration
		nConfigurations[s][(cal,car,cbl,cbr)] += nExpectedConfigurations - configurationsReported

	if ("configurations" in debug):
		dbg_report_configurations(nConfigurations)

	return nConfigurations


def configurations_count(hashA,hashB,windowSize):
	w = windowSize
	return (len(hashA)-w) * (len(hashB)-w)


unreportedConfiguration = (0,2,2,2,2)

def configurations(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize):
	# faster than all_configurations(), but harder to prove correct
	w = windowSize
	iMax = len(hashA)-(w+1)
	jMax = len(hashB)-(w+1)

	hasBeenReported = set()
	for ii in aPosToBPos:
		jj = aPosToBPos[ii]
		for i in range(ii-w,ii+1):           # min is ii-w, max is ii
			if (i < 0) or (i > iMax): continue
			for j in range(jj-w,jj+1):       # min is aPosToBPos[i]-w, max is aPosToBPos[i]
				if (j < 0) or (j > jMax): continue
				if ((i,j) in hasBeenReported): continue
				s = shared[(i+1,j+1)] if ((i+1,j+1) in shared) else 0
				(cal,car,cbl,cbr) = configuration_matrix(hashA,hashB,aPosToBPos,bPosToAPos,w,i,j)
				yield (s,cal,car,cbl,cbr,i,j)
				hasBeenReported.add((i,j))


def all_configurations(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize):
	# easier to prove correct than configurations(), but slower
	w = windowSize
	for i in range(len(hashA)-w):        # min is 0, max is |a|-(w+1)
		for j in range(len(hashB)-w):    # min is 0, max is |b|-(w+1)
			s = shared[(i+1,j+1)] if ((i+1,j+1) in shared) else 0
			(cal,car,cbl,cbr) = configuration_matrix(hashA,hashB,aPosToBPos,bPosToAPos,w,i,j)
			yield (s,cal,car,cbl,cbr,i,j)


def configuration_matrix(hashA,hashB,aPosToBPos,bPosToAPos,windowSize,i,j):
	w = windowSize

#	if ("configsnoop" in debug) and (i == 2) and (j == 3):
#		print("~~~~~",file=stderr)
#		print("hashA[%d] == hashB[%d] is %s" % (i,j,hashA[i] == hashB[j]),file=stderr)
#		if (i   not in aPosToBPos): print("match for A[%d] : none" % i,file=stderr)
#		else:                       print("match for A[%d] : %d <= %d <= %d is %s" % (i,j+1,aPosToBPos[i],j+w,j+1 <= aPosToBPos[i] <= j+w),file=stderr)
#		print("hashA[%d] == hashB[%d] is %s" % (i+w,j+w,hashA[i+w] == hashB[j+w]),file=stderr)
#		if (i+w not in aPosToBPos): print("match for A[%d] : none" % (i+w),file=stderr)
#		else:                       print("match for A[%d] : %d <= %d <= %d is %s" % (i+w,j+1,aPosToBPos[i+w],j+w-1,j+1 <= aPosToBPos[i+w] <= j+w-1),file=stderr)
#		print("hashB[%d] == hashA[%d] is %s" % (j,i,hashB[j] == hashA[i]),file=stderr)
#		if (j   not in bPosToAPos): print("match for B[%d] : none" % j,file=stderr)
#		else:                       print("match for B[%d] : %d <= %d <= %d is %s" % (j,i+1,bPosToAPos[j],i+w,i+1 <= bPosToAPos[j] <= i+w),file=stderr)
#		print("hashB[%d] == hashA[%d] is %s" % (j+w,i+w,hashB[j+w] == hashA[i+w]),file=stderr)
#		if (j+w not in bPosToAPos): print("match for B[%d] : none" % (j+w),file=stderr)
#		else:                       print("match for B[%d] : %d <= %d <= %d is %s" % (j+w,i+1,bPosToAPos[j+w],i+w-1,i+1 <= bPosToAPos[j+w] <= i+w-1),file=stderr)
#		print("~~~~~",file=stderr)

	if   (hashA[i]   == hashB[j]):          cal = 0
	elif (i not in aPosToBPos):             cal = 2
	elif (j+1 <= aPosToBPos[i]   <= j+w):   cal = 1
	else:                                   cal = 2

	if   (hashA[i+w] == hashB[j+w]):        car = 0
	elif (i+w not in aPosToBPos):           car = 2
	elif (j+1 <= aPosToBPos[i+w] <= j+w-1): car = 1
	else:                                   car = 2

	if   (hashB[j]   == hashA[i]):          cbl = 0
	elif (j not in bPosToAPos):             cbl = 2
	elif (i+1 <= bPosToAPos[j]   <= i+w):   cbl = 1
	else:                                   cbl = 2

	if   (hashB[j+w] == hashA[i+w]):        cbr = 0
	elif (j+w not in bPosToAPos):           cbr = 2
	elif (i+1 <= bPosToAPos[j+w] <= i+w-1): cbr = 1
	else:                                   cbr = 2

	return (cal,car,cbl,cbr)


# hash_sequence--
#	Convert a sequence to hashed kmers

ntToBits = {"A":0 , "C":1, "G":2, "T":3}

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

		kmerBits = 0
		for nt in kmer:
			kmerBits = (kmerBits<<2) + ntToBits[nt]

		hashSeq += [hashFunc(kmerBits)]

	return tuple(hashSeq)


# shared_hash_positions--
#	Return a dict that maps positions in hash sequence A to positions in B that
#	have the same hash.

def shared_hash_positions(hashA,hashB,duplicates=None):
	hashToB = {}
	for (bPos,hB) in enumerate(hashB):
		if (hB in hashToB):
			# sequences are not supposed to have any duplciates
			if (duplicates == None): raise ValueError  
			duplicates += [("B",hB,hashToB[hB],bPos)]
		hashToB[hB] = bPos

	if (duplicates != None) and (duplicates != []):
		raise ValueError

	aPosToBPos = {}
	hashToA = {}
	for (aPos,hA) in enumerate(hashA):
		if (hA in hashToA):
			# sequences are not supposed to have any duplciates
			if (duplicates == None): raise ValueError  
			duplicates += [("A",hA,hashToA[hA],aPos)]
		hashToA[hA] = aPos
		if (hA in hashToB): aPosToBPos[aPos] = hashToB[hA]

	if (duplicates != None) and (duplicates != []):
		raise ValueError

	return aPosToBPos


# jaccard_by_kmers--
#	Compute the true jaccard between two sequences, from their kmer sets

class KmerData: pass

def jaccard_by_kmers(seqA,seqB,kmerSize,canonical=False):
	kmersA = distinct_kmers(seqA,kmerSize,canonical=canonical)
	kmersB = distinct_kmers(seqB,kmerSize,canonical=canonical)

	numKmersA = len(kmersA)
	numKmersB = len(kmersB)

	nIntersection = len(kmersA.intersection(kmersB))
	nUnion = numKmersA + numKmersB - nIntersection

	if ("kmers" in debug):
		dbg_report_kmers(kmersA,kmersB)

	kd = KmerData()
	kd.jaccard       = nIntersection / nUnion   # (non-truncating divide)
	kd.nIntersection = nIntersection
	kd.nUnion        = nUnion
	return kd


# distinct_kmers--
#	return the set of distinct kmers in a sequence

def distinct_kmers(seq,kmerSize,canonical=False):
	seqLen = len(seq)
	distinct = set()

	if (canonical):
		revSeq = reverse_complement(seq)

	for ix in range(kmerSize,seqLen+1):
		kmer = seq[ix-kmerSize:ix]
		if (canonical):
			revKmer = revSeq[seqLen-ix:seqLen-ix+kmerSize]
			kmer = min(kmer,revKmer)
		distinct.add(kmer)

	return distinct


# jaccard_by_minimizers--
#	Compute the jaccard estimate between two hash sequences, from their
#	minimizer sets

class MinimizerData: pass

def jaccard_by_minimizers(hashA,hashB,windowSize,winnower):
	minimizersA = set([v for (v,ix) in winnower(hashA,windowSize)])
	minimizersB = set([v for (v,ix) in winnower(hashB,windowSize)])

	nMinimizersA = len(minimizersA)
	nMinimizersB = len(minimizersB)

	nIntersection = len(minimizersA.intersection(minimizersB))
	nUnion = nMinimizersA + nMinimizersB - nIntersection

	md = MinimizerData()
	if (nUnion == 0): md.jaccard = 0.0
	md.jaccard        = 0.0 if (nUnion == 0) \
	                    else (nIntersection / nUnion)  # (non-truncating divide)
	md.nIntersection  = nIntersection
	md.nUnion         = nUnion
	md.nMinimizersA   = nMinimizersA
	md.nMinimizersB   = nMinimizersB
	return md


# fasta_pairs--

def fasta_pairs(f1,f2=None,f2Name=None):
	f2Generator = fasta_sequences(f2) if (f2 != None) else None

	seqA = seqB = None
	for (seqName,seq) in fasta_sequences(f1):
		if (f2Generator == None):
			# one input file; alternate sequences from f1 as seqA and seqB
			if (seqA == None):
				(nameA,seqA) = (seqName,seq.upper())
				continue
			(nameB,seqB) = (seqName,seq.upper())
		else:
			# two input files; seqA is from f1 and seqB from f2
			(nameA,seqA) = (seqName,seq.upper())
			try:
				(nameB,seqB) = next(f2Generator)
			except StopIteration:
				print("WARNING: \"%s\" doesn't contain enough sequences" % f2Name,file=stderr)
				print("%s and any sequences following it had no sequence to compare to" % nameA,file=stderr)
				break
			seqB = seqB.upper()

		yield (nameA,seqA,nameB,seqB)

		seqA = seqB = None

	if (f2Generator == None):
		# one input file; make sure it had an even number of sequences
		if (seqA != None):
			print("WARNING: input contains an odd number of sequences",file=stderr)
			print("%s had no sequence to compare to" % nameA,file=stderr)
	else:
		# two input files; make sure f2 has no leftover sequences
		try:
			(nameB,seqB) = next(f2Generator)
			print("WARNING: \"%s\" contains too many sequences" % f2Name,file=stderr)
			print("%s and any sequences following it had no sequence to compare to" % nameB,file=stderr)
		except StopIteration:
			pass


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


# debug prints

def dbg_report_kmers(kmersA,kmersB):
	print("A: {%s}" % ",".join(kmersA),file=stderr)
	print("B: {%s}" % ",".join(kmersB),file=stderr)

def dbg_report_kmer_seqs(seqA,seqB,kmerSize,canonical=False):
	print("A: {%s}" % ",".join(["%d:%s"%(ix,kmer)for (ix,kmer) in enumerate(dbg_kmerize(seqA,kmerSize,canonical=canonical))]),file=stderr)
	print("B: {%s}" % ",".join(["%d:%s"%(ix,kmer) for (ix,kmer) in enumerate(dbg_kmerize(seqB,kmerSize,canonical=canonical))]),file=stderr)

def dbg_kmerize(seq,kmerSize,canonical=False):
	for ix in range(kmerSize,len(seq)+1):
		kmer = seq[ix-kmerSize:ix]
		if (canonical):
			revKmer = reverse_complement(kmer)
			kmer = min(kmer,revKmer)
		yield kmer

def dbg_report_hashes(hashA,hashB,kmerSize):
	nybbles = (2*kmerSize+3) // 4
	print("A: {%s}" % ",".join(["%d:%0*X"%(ix,nybbles,v) for (ix,v) in enumerate(hashA)]),file=stderr)
	print("B: {%s}" % ",".join(["%d:%0*X"%(ix,nybbles,v) for (ix,v) in enumerate(hashB)]),file=stderr)

def dbg_report_shared(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize):
	w = windowSize
	print("A: {%s}" % ",".join(["%d->%d"%(ix,aPosToBPos[ix]) for ix in aPosToBPos]),file=stderr)
	print("B: {%s}" % ",".join(["%d->%d"%(ix,bPosToAPos[ix]) for ix in bPosToAPos]),file=stderr)
	keysConsidered = set()
	for i in range(len(hashA)-(w-1)):        # min is 0, max is |a|-w
		for j in range(len(hashB)-(w-1)):
			keysConsidered.add((i,j))
			if ((i,j) not in shared): continue
			print("S(%d,%d,10) = %d" % (i,j,shared[(i,j)]),file=stderr)
	for (i,j) in shared:
		if ((i,j) not in keysConsidered):
			print("S(%d,%d,10) = %d (problem?)" % (i,j,shared[(i,j)]),file=stderr)

def dbg_report_configuration(i,j,s,cal,car,cbl,cbr):
	print("config(%d,%d) = [%d,%d;%d,%d;%d]"
	    % (i,j,cal,car,cbl,cbr,s),
	      file=stderr)

def dbg_report_configurations(nConfigurations):
	sValues = list(nConfigurations.keys())
	sValues.sort()
	sValues.reverse()
	for s in sValues:
		cValues = list(nConfigurations[s].keys())
		cValues.sort()
		for (cal,car,cbl,cbr) in cValues:
			print("N([%d,%d;%d,%d;%d]) --> %d"
			    % (cal,car,cbl,cbr,s,
			       nConfigurations[s][(cal,car,cbl,cbr)]),
			      file=stderr)

def dbg_report_alpha(alpha,s):
	for t in [0,1,2]:
		print("alpha_%d(%d) = %d" % (t,s,alpha[t]),file=stderr)

def dbg_report_CI(dbgCI):
	for (s,deltaCI,parts) in dbgCI:
		prefix = "CI =" if (s == 0) else "+"
		partsStr = "" if (len(parts) == 0) \
		           else (", %s" % ("+".join(["%0.6f*%s"%(part[0],part[1]) for part in parts])))
		print("%s %0.6f (for s=%d%s)" % (prefix,deltaCI,s,partsStr),file=stderr)


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
