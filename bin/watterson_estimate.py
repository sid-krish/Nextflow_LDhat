#!/usr/bin/env python
from math import log
import pandas as pd
import sys


def watterson_estimate(segregating_sites, genome_len, samples):
    # Theta per site
    k = 1
    calc_sum = 0
    while k < samples:
        calc_sum += 1 / k
        k += 1

    sum_inverse = calc_sum ** -1

    calc_log = log(genome_len / (genome_len - segregating_sites))

    theta = sum_inverse * calc_log

    return theta


# def output_theta(theta):
#     with open('theta.txt', 'w') as file:
#         file.write(f"{theta}")


if __name__ == '__main__':
    num_segregating_sites = int(sys.argv[1])
    genome_len = int(sys.argv[2])
    sample_size = int(sys.argv[3])


    theta = watterson_estimate(num_segregating_sites, genome_len, sample_size)
    # output_theta(theta)

    sys.stdout.write(str(round(theta,5))) # for next steps
