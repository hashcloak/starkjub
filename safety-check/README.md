Starkjub supporting evidence
--------------------------

Disclaimer: This repository is forked from [daira/jubjub](https://github.com/daira/jubjub), and adapted to new twisted Edwards curve, Starkjub, thanks to ZCash and Daira to make [daira/jubjub](https://github.com/daira/jubjub) open-source

This repository contains supporting evidence that the twisted Edwards curve
ax^2 + y^2 = 1 + dx^2.y^2, where a = 146640 and d = 146636
GF= (3618502788666131213697322783095070105623107215331596699973092056135872020481), which is the prime field used in [STARK Curve](https://docs.starkware.co/starkex/crypto/stark-curve.html)
satisfies the [SafeCurves criteria](https://safecurves.cr.yp.to/index.html).

The script ``verify.sage`` is based on
[this script from the SafeCurves site](https://safecurves.cr.yp.to/verify.html),
modified

* to support twisted Edwards curves;
* to generate a file 'primes' containing the primes needed for primality proofs,
  if it is not already present;
* to change the directory in which Pocklington proof files are generated
  (``proof/`` rather than ``../../../proof``), and to create that directory
  if it does not exist.

Prerequisites:

* apt-get install sagemath
* sage --pip install sortedcontainers

Run ``sage verify.sage .``, or ``./run.sh`` to also print out the results.

Note that the "rigidity" criterion cannot be checked automatically.
