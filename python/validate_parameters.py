#!/usr/bin/env python3
"""Fail before running expensive tools when settings are inconsistent."""
import os
import math


def main():
    integers = ['THREADS','MIN_DP','MAX_DP','LD_WINDOW','LD_STEP','N_PC','K_MIN',
                'K_MAX','N_REPS','CV_FOLDS','FST_WINDOW','FST_STEP','BLOCK_BP',
                'N_BOOT','MIN_CALLED','REPS','BLOCK_SIZE']
    for key in integers:
        value = int(os.environ[key])
        if value < 1:
            raise ValueError(key + ' must be a positive integer')
    for key, upper in [('MAF', .5), ('MISS', 1.), ('LD_R2', 1.)]:
        value = float(os.environ[key])
        if not math.isfinite(value) or not 0 < value <= upper:
            raise ValueError('Invalid ' + key)
    if not math.isfinite(float(os.environ['MIN_QUAL'])) or float(os.environ['MIN_QUAL']) < 0:
        raise ValueError('MIN_QUAL must be finite and nonnegative')
    for low,high in [('MIN_DP','MAX_DP'),('K_MIN','K_MAX'),('LD_STEP','LD_WINDOW'),('FST_STEP','FST_WINDOW')]:
        if int(os.environ[low]) > int(os.environ[high]):
            raise ValueError(low + ' must not exceed ' + high)
    if os.environ['APPLY_DP'] not in {'0','1'}:
        raise ValueError('APPLY_DP must be 0 or 1')
    if int(os.environ['MAX_M']) < 0 or not 0 <= int(os.environ['SEED']) < 2**31 - 1000000:
        raise ValueError('Invalid MAX_M or SEED')
    if int(os.environ['MIN_CALLED']) < 2 or int(os.environ['N_BOOT']) < 20:
        raise ValueError('MIN_CALLED >= 2 and N_BOOT >= 20 required')

if __name__ == '__main__':
    main()
