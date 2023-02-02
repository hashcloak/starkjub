# starkjub
A STARK-friendly elliptic curve designed to work within Cairo circuits on Starknet

## Overview
Starkjub is a Twisted edwards curve to be used within Cairo circuits for Starknet. This repository contains code from [Daira Hopwood](https://github.com/daira/jubjub) and [Pratyush Mishra](https://gist.github.com/Pratyush/68c17a5ed34bb4323ed1b4641f5c181c#file-sample_edwards-sage) in order to ensure that Starkjub defines a valid Twisted edwards curve that satisfies the SafeCurves criteria.

## Parameters
Starkjub is a TWisted edwards curve with a = 146640 and d = 146636 over the finite field GF(3618502788666131213697322783095070105623107215331596699973092056135872020481). 

The values for a and d were found from the procedure described [here](https://www.rfc-editor.org/rfc/rfc7748#appendix-A.1) and converted from the Montgomery form. The prime used is the Stark prime, defined [here](https://docs.starkware.co/starkex/crypto/stark-curve.html)

## Install Instructions
Make sure you have [sage](https://www.sagemath.org/) and Python installed on your system.

```
git clone git@github.com:hashcloak/starkjub.git
cd src
python starkjub.py
```
If you want to run the safecurves check yourself then, `cd safety-check` and follow the instructions in the README.
