#!/usr/bin/env python3
"""
Compute the Jaccard correction for pairs of sequences, relating to winnowed
minimizers.

This is a snapshot of the program, with internal formulas relating to an
earlier version of the associated manuscript.
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
                                splitmix64_hash_names
from winnowed_minimizers import winnowed_minimizers_linear

programName    = "jaccard_correction_test"
programVersion = "0.4.0"


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
                          type is either minimap2, murmurhash3, or splitmix64
                          (default is minimap2.0)
  --prng=<string>         set seed for PRNG; NOT to be confused with the hash
                          seed; this may include substrings "{date}" and
                          "{time}"
  --report:configs        include configuration counts in the output; this
                          greatly increases the number of output columns
  --inhibit:correction    don't compute any of the statistics related to the
                          jaccard correction
  --head=<number>         limit the number of input sequence pairs
  --progress=<number>     periodically report how many sequence pairs we've
                          processed
  --version               report this program's version number

Compute the Jaccard correction for pairs of sequences, relating to winnowed
minimizers.

NOTE: This is a snapshot of the program, with internal formulas relating to an
      earlier version of the associated manuscript.

If fasta_file_2 is provided, the first sequence in fasta_file_1 is compared to
the first from fasta_file_2, the second is compared to the second, and so on.

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
	reportConfigs     = False
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
		elif (arg.startswith("--prng=")):
			if (argVal.lower() == "none"):
				prngSeed = None
			else:
				prngSeed = argVal
				prngSeed = prngSeed.replace("{date}",strftime("%d/%m/%Y"))
				prngSeed = prngSeed.replace("{time}",strftime("%I/%M/%S"))
		elif (arg in ["--report:configs","--report:configurations"]):
			reportConfigs = True
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
	#
	# variable                              name in header      name in manuscript
	# --------                              --------------      ------------------
	# nameA                                 nameA            
	# nameB                                 nameB            
	# replicateNum                          rep
	# hashSeed                              hash.seed
	# numReplicates                         replicates
	# windowSize                            w                   w
	# kmerSize                              k                   k
	# len(seqA)                             length.nt           (L+k-1)
	# len(hashA)                            |a|                 L
	# kd.nIntersection                      I(A,B)              I(A,B)
	# kd.nUnion                             U(A,B)              U(A,B) 
	# kd.jaccard                            J(A,B)              J(A,B)
	# md.nIntersection or nIntersectionAvg  I(A,B;w)            Ihat(A,B;w) or average(Ihat(A,B;w))
	# md.nUnion or nUnionAvg                U(A,B;w)            Uhat(A,B;w) or average(Uhat(A,B;w))
	# md.jaccard or jaccardAvg              J(A,B;w)            Jhat(A,B;w) or Jbar(A,B;w)
	# cd.scriptD                            D(A,B;w)            script D(A,B;w)
	# cd.jaccardFromD                       Jd(A,B;w)
	# cd.scriptC                            C(A,B;w)            script C(A,B;w)
	# cd.bias                               Bias(A,B;w)         script B(A,B;w)

	showAllReplicates = ("replicates" in debug)

	header =  ["nameA","nameB"]
	if (numReplicates > 1):
		if (showAllReplicates): header += ["rep","hash.seed"]
		else:                   header += ["replicates"]
	header += ["w","k","length.nt","|a|"]
	header += ["I(A,B)","U(A,B)","J(A,B)"]
	if ("density" in debug): header += ["2L/(w+1)","|A.min(w)|","|B.min(w)|"]
	header += ["I(A,B;w)","U(A,B;w)","J(A,B;w)"]
	header += ["D(A,B;w)","Jd(A,B;w)","C(A,B;w)","Bias(A,B;w)"]
	header += ["J(A,B;w)-J(A,B)","I(A,B;w)-C(A,B;w)"]
	if (reportConfigs):
		#viableConfigurations = []
		#for (cal,car,cbl,cbr) in cartesian_product([0,1,2],[0,1,2],[0,1,2],[0,1,2]):
		#	viableConfigurations += [(cal,car,cbl,cbr)]
		viableConfigurations = [(0,0,0,0),(0,1,0,1),(0,1,0,2),(0,2,0,1),(0,2,0,2),
                                (2,0,2,0),(2,1,2,1),(2,1,2,2),(2,2,2,1),(2,2,2,2),
                                (2,1,1,1),(2,2,1,1),(1,1,2,1),(1,1,2,2),(1,0,1,0),
                                (1,0,2,0),(2,0,1,0)]
		for (cal,car,cbl,cbr) in viableConfigurations:
			for s in range(windowSize+1):             # (min is 0, max is w)
				header += ["N(%d,%d;%d,%d;%d)"%(cal,car,cbl,cbr,s)]

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
		# kd.nIntersection is manuscript's Sp^k(A) intersect Sp^k(B)
		# kd.nUnion        is manuscript's Sp^k(A) union Sp^k(B)
		# kd.jaccard       is manuscript's J(A,B)

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
			# md.nIntersection is manuscript's I(A,B;w)
			# md.nUnion        is manuscript's U(A,B;w)
			# md.jaccard       is manuscript's J(A,B;w)

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
				# cd.scriptD is manuscript's script D(A,B;w)
				# cd.scriptC is manuscript's script C(A,B;w)
				# cd.bias    is manuscript's script B(A,B;w)

				if ("C" in debug):
					dbg_report_script_c(cd.dbgScriptC)

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
				line += ["%d\t%.6f\t%.6f\t%.6f" % (cd.scriptD,cd.jaccardFromD,cd.scriptC,cd.bias)]
			else:
				line += ["NA\tNA\tNA\tNA"]
			if (numReplicates == 1) or (showAllReplicates):
				if (computeCorrection):
					line += ["%.6f\t%.6f"   % (md.jaccard-kd.jaccard,md.nIntersection-cd.scriptC)]
				else:
					line += ["%.6f\tNA"     % (md.jaccard-kd.jaccard)]
			elif (replicateNum == numReplicates-1):
				if (computeCorrection):
					line += ["%.6f\t%.6f"   % (jaccardAvg-kd.jaccard,nIntersectionAvg-cd.scriptC)]
				else:
					line += ["%.6f\tNA"     % (jaccardAvg-kd.jaccard)]

			if (reportConfigs):
				for (cal,car,cbl,cbr) in viableConfigurations:
					for s in range(windowSize+1):             # (min is 0, max is w)
						if (s in cd.nConfigurations) and ((cal,car,cbl,cbr) in cd.nConfigurations[s]):
							n = cd.nConfigurations[s][(cal,car,cbl,cbr)]
						else:
							n = 0
						line += [str(n)]

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
	# shared[(i,j)] is manuscript's S(i,j,w), but only when S(i,j,w)>0
	# nConfigurations[s][(cal,car,cbl,cbr)] is manuscript's N[(cal,car;cbl,cbr;s)]

	shared = window_shared_counts(hashA,hashB,windowSize,aPosToBPos)

	if ("shared" in debug):
		dbg_report_shared(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize)

	nConfigurations = configuration_counts(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize)

	L = len(hashA)
	w = windowSize

	if ("C" in debug):
		dbgScriptC = []

	scriptD = 0
	for s in range(windowSize+1):
		if (s not in nConfigurations): continue
		for (cal,cbl) in cartesian_product([0,1,2],[0,1,2]):
			scriptD += nConfigurations[s][(cal,0,cbl,0)]
	jaccardFromD = float(scriptD) / (2*L - scriptD)

	scriptC = 0                                       # scriptC is manuscript's script C(A,B;w)
	for s in range(windowSize+1):                     # (min is 0, max is w)
		nOmega = {}                                   # nOmega[t] is manuscript's N(omega_{t,s})
		for t in [0,2]:                               # (note that N(omega_{1,s}) is not needed)
			nOmega[t] = 0
			for (car,cbr) in [(1,2),(2,1),(0,0)]:
				nOmega[t] += nConfigurations[s][(t,car,t,cbr)]
		if ("omega" in debug):
			dbg_report_n_omega(nOmega,s)

		g = g_func                                    # g_func(w,s,alpha,beta) is manuscript's g(s,alpha,beta)
		scriptC +=     nConfigurations[s][(0,2,0,2)] * g(w,s,0,1)
		scriptC += 2 * nConfigurations[s][(2,2,2,2)] * g(w,s,0,2)
		scriptC += 2 * nConfigurations[s][(2,1,2,1)] * g(w,s,2,2)
		for (cal,cbl) in cartesian_product([0,1,2],[0,1,2]):
			scriptC += nConfigurations[s][(cal,0,cbl,0)] * g(w,s,s-1,0)
		for (cal,car,cbl) in cartesian_product([0,2],[0,1,2],[0,1,2]):
			scriptC += nConfigurations[s][(cal,car,cbl,1)] * g(w,s,s-1,1)
		for (cal,cbl,cbr) in cartesian_product([0,1,2],[0,2],[0,1,2]):
			scriptC += nConfigurations[s][(cal,1,cbl,cbr)] * g(w,s,s-1,1)
		if ("0.3.0" not in debug):
			# this term was left out of an earlier version of the manuscript,
			# and this program
			scriptC += nConfigurations[s][(0,1,0,1)] * g(w,s,2,1)
		scriptC +=     nOmega[0] * g(w,s,1,1)
		scriptC += 2 * nOmega[2] * g(w,s,1,2)

		if ("C" in debug):
			deltaScriptC = scriptC if (s==0) else (scriptC-prevScriptC)
			prevScriptC = scriptC
			parts = []
			parts += [(  nConfigurations[s][(0,2,0,2)],g(w,s,0,1))]
			parts += [(2*nConfigurations[s][(2,2,2,2)],g(w,s,0,2))]
			parts += [(2*nConfigurations[s][(2,1,2,1)],g(w,s,2,2))]
			parts += [(sum([nConfigurations[s][(cal,0,cbl,0)] for (cal,cbl) in cartesian_product([0,1,2],[0,1,2])]),g(w,s,s-1,0))]
			parts += [(sum([nConfigurations[s][(cal,car,cbl,1)] for (cal,car,cbl) in cartesian_product([0,2],[0,1,2],[0,1,2])]),g(w,s,s-1,1))]
			parts += [(sum([nConfigurations[s][(cal,1,cbl,cbr)] for (cal,cbl,cbr) in cartesian_product([0,1,2],[0,2],[0,1,2])]),g(w,s,s-1,1))]
			if ("0.3.0" not in debug):
				parts += [(nConfigurations[s][(0,1,0,1)],g(w,s,2,1))]
			parts += [(  nOmega[0],g(w,s,1,1))]
			parts += [(2*nOmega[2],g(w,s,1,2))]
			dbgScriptC += [(s,deltaScriptC,parts)]

	bias = (scriptC/((float(4*L)/(w+1))-scriptC)) - (float(scriptD)/(2*L-scriptD))
	# bias is manuscript's script B Bias(A,B;w)

	cd = CorrectionData()
	cd.scriptD         = scriptD
	cd.jaccardFromD    = jaccardFromD
	cd.scriptC         = scriptC
	cd.bias            = bias
	cd.nConfigurations = nConfigurations
	if ("C" in debug): cd.dbgScriptC = dbgScriptC
	return cd


# window_shared_counts
#	return shared s.t. shared[(i,j)] is the count of shared hash values in
#	A[i..i+w-1] and B[j..j+w-1]. shared is only populated for (i,j) that lead
#	to non-zero counts.

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
# The vast majority of configs are [2,2;2,2;0]. configurations() avoids
# reporting most of those, but computes all the others. The count for the
# unreported [2,2;2,2;0] is then deduced from the difference in the number of
# configurations reported and the expected number, configurations_count().
# Note that configurations() may report *some* [2,2;2,2;0].
#
# Hypothesis:
#	All non-[2,2;2,2;0] configurations fall into this set:
#	union over i',j' s.t. a[i']=b[j'] of {i,j s.t. i<=i'-w<=i and j<=j'-w<=j}
#
# The config at (i,j) is *not* [2,2;2,2;0] if any of these are true
#	s!=0     S(i+1,j+1,w) != 0           (i+1,j+1) in shared
#	cal!=2   a[i] is in b[j..j+w]        j   <= aPosToBPos[i]   <= j+w
#	car!=2   a[i+w] is in b[j+1..j+w]    j+1 <= aPosToBPos[i+w] <= j+w
#	cbl!=2   b[j] is in a[i..i+w]        i   <= bPosToAPos[j]   <= i+w
#	cbr!=2   b[j+w] is in a[i+1..i+w]    i+1 <= bPosToAPos[j+w] <= i+w

def configuration_counts(hashA,hashB,aPosToBPos,bPosToAPos,shared,windowSize):
	configurations_reporter = all_configurations if ("allconfigs" in debug) \
	                     else configurations

	nExpectedConfigurations = configurations_count(hashA,hashB,windowSize)

	# nConfigurations[s][(cal,car,cbl,cbr)] is manuscript's N[(cal,car;cbl,cbr;s)]
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


# g_func--
#	This is g(s,alpha,beta) from the manuscript.
#
# nota bene: We back this function with a cache. This implementation would
#            probably be cleaner if it used the python memoization paradigm.

g_func_hash = {}
def g_func(w,s,alpha,beta,blindToCache=False):
	cacheKey = (w,s,alpha,beta)
	if (not blindToCache) and (cacheKey in g_func_hash):
		return g_func_hash[cacheKey]

	g = s - alpha
	for i in range(0,beta+1):    # min is 0, max is beta
		g /= 2*w-s+i

	g_func_hash[cacheKey] = g
	return g


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
	assert (len(seqA) > 0) or (len(seqB) > 0)

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
			print("S(%d,%d,%d) = %d" % (i,j,windowSize,shared[(i,j)]),file=stderr)
	for (i,j) in shared:
		if ((i,j) not in keysConsidered):
			print("S(%d,%d,%d) = %d (problem?)" % (i,j,windowSize,shared[(i,j)]),file=stderr)

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

def dbg_report_n_omega(nOmega,s):
	for t in [0,2]:
		print("N(omega_{%d,%d} = %d" % (t,s,nOmega[t]),file=stderr)

def dbg_report_script_c(dbgScriptC):
	for (s,deltaScriptC,parts) in dbgScriptC:
		prefix = "scriptC =" if (s == 0) else "+"
		partsStr = "" if (len(parts) == 0) \
		           else (", %s" % ("+".join(["%s*%0.6f"%(part[0],part[1]) for part in parts])))
		print("%s %0.6f (for s=%d%s)" % (prefix,deltaScriptC,s,partsStr),file=stderr)


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
