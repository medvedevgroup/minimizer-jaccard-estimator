#!/usr/bin/env python3

from collections import deque

# winnowed_minimizers--
#	Yields a sequence of (v,ix) pairs; v is a minimizer value, ix is it's
#	location in perm.
#
# note: we expect values in perm to be unique

def winnowed_minimizers(perm,windowSize,circular=False):
	if (circular):
		for mini in winnowed_minimizers_circular(perm,windowSize): yield mini
	else:
		for mini in winnowed_minimizers_linear(perm,windowSize): yield mini


def winnowed_minimizers_linear(perm,windowSize):
	# note: we expect values in perm to be unique
	# note: we treat perm as a linear sequence
	# note: if len(perm) < windowSize, no minimizers will be reported

	history = deque()	# ordered oldest (front) to most recent (back)
						# this will contain (value,position) pairs

	minimizers = set()
	for ix in range(len(perm)):
		windowIx = ix - (windowSize-1)
		v = perm[ix]

		# (at the back) remove anything worse than the current value

		while (len(history) > 0) and (history[-1][0] > v):
			history.pop()

		# push the current value and its position to back of the queue

		history.append((v,ix))

		# (at the front) remove any minima that are no longer in the current
		# window

		while (len(history) > 0) and (history[0][1] < windowIx):
			history.popleft()

		# copy the minimizer from window, but only if we are seeing it for
		# the first time

		if (windowIx >= 0) and (len(history) > 0):
			vHistory = history[0]
			if (vHistory not in minimizers):
				yield vHistory
				minimizers.add(vHistory)


def winnowed_minimizers_circular(perm,windowSize):
	# note: we expect values in perm to be unique
	# note: we treat perm as a circular sequence

	history = deque()	# ordered oldest (front) to most recent (back)
						# this will contain (value,position) pairs

	minimizers = set()
	for ix in range(len(perm)+windowSize-1):
		windowIx = ix - (windowSize-1)
		wrapround = (ix >= len(perm))
		v = perm[ix] if (not wrapround) else perm[ix-len(perm)]

		# (at the back) remove anything worse than the current value

		while (len(history) > 0) and (history[-1][0] > v):
			history.pop()

		# push the current value and its position to back of the queue

		history.append((v,ix))

		# (at the front) remove any minima that are no longer in the current
		# window

		while (len(history) > 0) and (history[0][1] < windowIx):
			history.popleft()

		# copy the minimizer from window, but only if we are seeing it for
		# the first time

		if (windowIx >= 0) and (len(history) > 0):
			vHistory = history[0]
			seenBefore = (vHistory in minimizers)
			if (not seenBefore) and (vHistory[1] >= len(perm)):
				vHistory2 = (vHistory[0],vHistory[1]-len(perm))
				if (vHistory2 in minimizers):
					seenBefore = True
			if (not seenBefore):
				yield vHistory
				minimizers.add(vHistory)

