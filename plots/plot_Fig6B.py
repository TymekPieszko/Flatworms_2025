from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

# Example command
# python plot_Fig6B.py /data/biol-bdelloids/scro4331/Flatworms_2025/simulations/sim_output/fitness /data/biol-bdelloids/scro4331/Flatworms_2025/plots/Fig6B.png
in_dir = Path(sys.argv[1])
plot_file = sys.argv[2]


def get_gc_rate(dir):
    return float(dir.name.split("~")[-1])


asex_dirs = [d for d in in_dir.glob(f"SEX~0.0/GC~*") if d.is_dir()]
asex_dirs = sorted(asex_dirs, key=get_gc_rate)
sex_dirs = [d for d in in_dir.glob(f"SEX~1.0/GC~*") if d.is_dir()]
sex_dirs = sorted(sex_dirs, key=get_gc_rate)
print(asex_dirs)
### Calculate fitness under different asex scenarios
asex_values = []
for rate_dir in asex_dirs:
    reps = []
    for file in rate_dir.glob("*.txt"):
        values = np.loadtxt(file)
        print(values)
        if values.shape != (50,):
            continue
        reps.append(values[-1])
    asex_values.append((np.mean(reps), 1.96 * np.std(reps) / np.sqrt(len(reps))))

### Calculate fitness for the sexual scenario
sex_reps = []
for file in in_dir.rglob("SEX~1.0/GC~*/*.txt"):
    values = np.loadtxt(file)
    if values.shape != (50,):
        continue
    sex_reps.append(values[-1])

sex_mean = np.mean(sex_reps)
# sex_ci = 1.96 * np.std(sex_reps) / np.sqrt(len(sex_reps))

fix, ax = plt.subplots(figsize=(19, 13))
xpos = np.arange(len(asex_values))
plt.axhline(
    y=sex_mean,
    color="dimgray",
    linewidth=4,
    linestyle="--",
)
# plt.axhline(
#     y=sex_mean - sex_ci,
#     color="black",
#     linewidth=2,
#     linestyle="--",
#     alpha=0.5,
# )
# plt.axhline(
#     y=sex_mean + sex_ci,
#     color="black",
#     linewidth=2,
#     linestyle="--",
#     alpha=0.5,
# )
colors = [plt.cm.plasma(i / len(asex_values)) for i in range(len(asex_values))]
for i, color in enumerate(colors):
    plt.errorbar(
        x=xpos[i],
        y=asex_values[i][0],
        yerr=asex_values[i][1],
        fmt="o",
        linewidth=4,
        capsize=8,
        markersize=10,
        color="black",
    )
plt.ylim(0.9, 1)
plt.xticks(
    xpos,
    [f"{(get_gc_rate(d)*5000):.{2}e}" for d in asex_dirs],
    rotation=26,
    size=31,
)
plt.yticks(size=31)
plt.ylabel("Mean fitness", size=37, labelpad=32)
plt.xlabel("GC rate", size=37, labelpad=17)
plt.tight_layout()
plt.savefig(plot_file)
