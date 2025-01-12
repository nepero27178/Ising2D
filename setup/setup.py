#!/usr/bin/python3

MIN_SIZE = 10
FIRST_EXCLUDED_SIZE = 160 # max_size is: first_excluded_size - step
STEP = 10
SIZES = range(MIN_SIZE, FIRST_EXCLUDED_SIZE, STEP)

# Maximum deviation from interval center
MAX_DEVIATION= 0.02

# TODO: make this file imported everywhere
