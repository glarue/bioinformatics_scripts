#!/usr/bin/env python3

import sys
import fileinput
import statistics as st

numbers = []

for l in fileinput.input():
    l = l.strip()
    try:
        numbers.append(float(l))
    except ValueError:
        continue

count = len(numbers)
max_v = max(numbers)
min_v = min(numbers)
sum_v = sum(numbers)
mean = st.mean(numbers)
median = st.median(numbers)
stdev = st.stdev(numbers)

properties = ['count', 'max', 'min', 'sum', 'mean', 'median', 'stdev']
values = [count, max_v, min_v, sum_v, mean, median, stdev]

for p, v in zip(properties, values):
    print(p + '\t' + str(v))

sys.exit(0)