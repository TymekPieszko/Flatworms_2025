from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Example command
# python plot_Fig6A.py /data/biol-bdelloids/scro4331/Flatworms_2025/plots/genotypes /data/biol-bdelloids/scro4331/Flatworms_2025/plots/Fig6A.png
in_dir = Path(sys.argv[1])
plot_file = sys.argv[2]


def get_gc_rate(dir):
    return float(dir.name.split("~")[-1])


# Collect directories for genotype files under asex scenarios
rate_dirs = [d for d in in_dir.glob(f"SEX~0.0/GC~*") if d.is_dir()]


fix, ax = plt.subplots(figsize=(19, 13))
colors = [plt.cm.plasma(i / len(rate_dirs)) for i in range(len(rate_dirs))]
handles = []
for i, rate_dir in enumerate(rate_dirs):
    gc_rate = get_gc_rate(rate_dir)
    print(f"GC rate: {gc_rate}")
    reps_del = []
    reps_neut = []
    for geno_file in rate_dir.glob("*.txt"):
        df = pd.read_csv(geno_file, sep="\t", header=None, names=["s", "hom"])
        df["s"] = -df[
            "s"
        ]  # Reversing the sign of s inverts the x axis without actually inverting it!

        ########################
        ### DELETERIOUS MUTS ###
        ########################
        df_del = df[df["s"] != 0.0]
        s_quants = np.quantile(df_del["s"], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        h = [
            np.mean(df_del["hom"][(df_del["s"] >= start) & (df_del["s"] < end)])
            for start, end in zip(s_quants[:-1], s_quants[1:])
        ]  # Homozygsoity between quantiles of s!
        reps_del.append(h)

        ####################
        ### NEUTRAL MUTS ###
        ####################
        df_neut = df[df["s"] == 0.0]
        h = np.mean(df_neut["hom"])  # A single homozygosity value for neutral mutations
        reps_neut.append(h)

    ### Plot deleterious
    reps_del = np.array(reps_del)
    mean_del = np.mean(reps_del, axis=0)
    ci_del = 1.96 * (np.std(reps_del, axis=0) / np.sqrt(len(reps_del)))
    # print(mean_across_reps)
    # print(ci)
    handle = plt.errorbar(
        np.linspace(0.5, 4.5, 5),
        mean_del,
        yerr=ci_del,
        linewidth=4,
        marker="o",
        capsize=8,
        markersize=10,
        color=colors[i],
    )
    handles.append(handle)

    ### Plot neutral
    reps_neut = np.array(reps_neut)
    mean_neut = np.mean(reps_neut)
    ci_neut = 1.96 * (np.std(reps_neut) / np.sqrt(len(reps_neut)))
    plt.errorbar(
        [0],
        [mean_neut],
        yerr=[ci_neut],
        linewidth=4,
        marker="o",
        capsize=8,
        markersize=10,
        color=colors[i],
    )
plt.legend(
    handles,
    [f"{(get_gc_rate(d)*5000):.2e}" for d in rate_dirs],
    loc="upper right",
    fontsize=22,
    reverse=True,
)
plt.ylim(0.0, 1.0)
plt.xlim(-0.5, 5.5)
plt.xticks(
    [0.0] + list(np.linspace(0.5, 4.5, 5)),
    ["Neutral", "0-20", "20-40", "40-60", "60-80", "80-100"],
    size=31,
)
# Plot the homozygosity of neutral mutations in S. mediterranea
plt.plot(
    0.0,
    0.315,
    "D",
    color="limegreen",
    markersize=14,
)
plt.yticks(size=31)
plt.ylabel("Homozygosity", size=37, labelpad=30)
plt.xlabel("Deleteriousness (%)", size=37, labelpad=17)
plt.tight_layout()
plt.savefig(plot_file)
