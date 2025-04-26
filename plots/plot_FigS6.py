from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

# Example command
# python plot_FigS6.py /data/biol-bdelloids/scro4331/Flatworms_2025/simulations/sim_output/fitness /data/biol-bdelloids/scro4331/Flatworms_2025/plots/FigS6.png
in_dir = Path(sys.argv[1])
plot_file = sys.argv[2]


def get_gc_rate(dir):
    return float(dir.name.split("~")[-1])


rate_dirs = [d for d in in_dir.glob(f"SEX~0.0/GC~*") if d.is_dir()]
rate_dirs = sorted(rate_dirs, key=get_gc_rate)


fig = plt.figure(figsize=(10, 10))
colors = [plt.cm.plasma(i / len(rate_dirs)) for i in range(len(rate_dirs))]
for i, dir in enumerate(rate_dirs):
    reps = []
    for file in dir.glob("*.txt"):
        values = np.loadtxt(file)
        if values.shape != (50,):
            continue
        reps.append(values)
    means = np.mean(reps, axis=0)
    ci = 1.96 * np.std(reps, axis=0) / np.sqrt(len(reps))
    x = list(range(0, 50000, 1000))
    plt.plot(x, means, "-o", linewidth=2, markersize=3, color=colors[i])
    plt.fill_between(x, means - ci, means + ci, color=colors[i], alpha=0.3)

plt.ylim(0.88, 1)
lines = plt.gca().get_lines()
plt.legend(
    lines[::-1],
    [f"{(get_gc_rate(d)*5000):.2e}" for d in rate_dirs][::-1],
    fontsize=18,
)
plt.xticks(list(range(0, 55000, 10000)), fontsize=24)
plt.yticks(fontsize=24)
plt.ylabel("Mean fitness", fontsize=30)
plt.xlabel("Time [generations]", fontsize=30)
plt.tight_layout()
plt.savefig(plot_file)
plt.yticks(fontsize=24)
plt.ylabel("Mean fitness", fontsize=26)
plt.xlabel("Time (generations)", fontsize=26)
plt.tight_layout()
plt.savefig(plot_file)
