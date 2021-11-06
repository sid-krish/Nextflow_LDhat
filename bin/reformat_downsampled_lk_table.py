#! /usr/bin/env python
import sys
import pandas as pd

if __name__ == '__main__':

    # Reformat downsampled lookup table to traditional format. Can then be used for LDhat

    downsampled_lk_table = sys.argv[1]
    num_samples = sys.argv[2]
    theta = sys.argv[3]
    ldpop_rho_range = sys.argv[4]

    # downsampled_lk_table = "../Lookup_tables/lk_downsampled_135.csv"
    # num_samples = '135'
    # theta = "0.01"
    # ldpop_rho_range = "101,100"

    steps, rho_range = ldpop_rho_range.split(',')

    df_downsampled = pd.read_csv(downsampled_lk_table)

    num_entries = df_downsampled.shape[0]
    entries = [i + 1 for i in range(num_entries)]

    df_downsampled["entry"] = entries
    df_downsampled.set_index('entry', inplace=True)

    df_downsampled[['00', '01', '10', '11']] = df_downsampled['00 01 10 11'].str.split(' ', expand=True)
    df_downsampled['#'] = '#'
    df_downsampled[':'] = ':'

    df_downsampled.drop('00 01 10 11', axis=1, inplace=True)

    final_cols = ['#', '00', '01', '10', '11', ':'] + df_downsampled.columns.to_list()[:-6]

    df_downsampled = df_downsampled[final_cols]

    with open("lookupTable.txt", 'w') as file_out:
        file_out.write(f"{num_samples} {num_entries}\n")
        file_out.write(f"{1} {theta}\n")
        file_out.write(f"{steps} {rho_range}\n")
        file_out.write("\n")
        file_out.write("\n")

    df_downsampled.to_csv('lookupTable.txt', mode='a', sep=' ', header=False)
