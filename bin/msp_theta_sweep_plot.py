import os

import numpy as np
import pandas as pd
import seaborn as sns

def collect_results_sweep_1(rho, theta, genome_size, sample_size, seed):
    # utilise numpy for combinations
    # https://www.kite.com/python/answers/how-to-get-all-element-combinations-of-two-numpy-arrays-in-python

    mesh_grid = np.array(np.meshgrid(rho, theta, genome_size, sample_size, seed))
    # reshape.(-1,2), -1 is unknown dimension, let numpy figure it out
    # z.reshape(-1, -1)
    # ValueError: can only specify one unknown dimension
    sweep_1_combinations = mesh_grid.T.reshape(-1, 5)

    # Load data into dataframe
    recom_est_results_dir = f"/Users/Sid/Documents/Github/run_ldhat/LDhat_Output_theta_sweep/"

    col_names = ["rho_sim", "theta_sim", "genome_size_sim", "sample_size_sim", "seed_sim", "theta_est"]

    df_recom_est_resutls = pd.DataFrame(columns=col_names)

    for rho, theta, genome_size, sample_size, seed in sweep_1_combinations:
        prepended_filename = f"rho_{rho}_theta_{theta}_genome_size_{int(genome_size)}_sample_size_{int(sample_size)}_seed_{int(seed)}_"
        with open(f"{recom_est_results_dir}{prepended_filename}msp_out_theta_est.csv", 'r') as results:
            results_unfiltered = results.readlines()[0].strip()

        to_append = [rho, theta, genome_size, sample_size, seed] + [results_unfiltered]
        df_recom_est_resutls.loc[len(df_recom_est_resutls)] = to_append

    return df_recom_est_resutls


if __name__ == '__main__':
    # Sweep 1: Recombination rate estimation
    rho_sweep_1 = [0.0] # unscaled r values. rho = 2 . p . N_e . r . tractlen
    theta_sweep_1 = [0.0005, 0.0025, 0.005]  # unscaled u values. theta = 2 . p . N_e . u
    genome_size_sweep_1 = [25000]
    sample_size_sweep_1 = [50,100,150,200]
    seed_sweep_1 = [1,2,3,4,5,6,7,8,9,10]

    recom_tract_len = 1000

    collected_results_sweep_1_df = collect_results_sweep_1(rho_sweep_1, theta_sweep_1, genome_size_sweep_1,
                                                           sample_size_sweep_1, seed_sweep_1)

    # process and export df for plotting
    # Since fastsimbac does 2n*r only scaling by tract len is needed

    collected_results_sweep_1_df["scaled_theta_sim"] = collected_results_sweep_1_df["theta_sim"].apply(
        lambda x: x * 2)

    # cols_to_drop = ["rho", "sample_size", "genome_size"]
    # collected_results_sweep_1_df.drop(columns=cols_to_drop, inplace=True)

    reorder_cols = ['rho_sim', 'theta_sim', 'scaled_theta_sim', 'genome_size_sim', 'sample_size_sim',
                    'seed_sim', "theta_est"]

    collected_results_sweep_1_df = collected_results_sweep_1_df.reindex(columns=reorder_cols)

    collected_results_sweep_1_df.to_csv("collected_results_msp_theta_sweep.csv", index=None)

    collected_results_sweep_1_df = collected_results_sweep_1_df.astype('float64')
    
    # Plot results
    sns.set_theme(style="whitegrid")

    ax = sns.boxplot(data=collected_results_sweep_1_df,x="scaled_theta_sim", y="theta_est", hue="sample_size_sim", palette="Greens")
    ax.legend(title='genomes')
    
    ax.set_title("LDhat")

    ax.set(ylabel="Estimated \u03B8", xlabel="Simulated \u03B8")

    ax.figure.savefig("collected_results_msp_theta_sweep.png",bbox_inches='tight', dpi=500)