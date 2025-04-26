import tskit, msprime, sys, tqdm
from pathlib import Path
import numpy as np

# Example command
# python calc_genotypes.py /data/biol-bdelloids/scro4331/Flatworms_2025/simulations/sim_output/ts/ /data/biol-bdelloids/scro4331/Flatworms_2025/plots/genotypes
in_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])


def get_gc_rate(dir):
    return float(dir.name.split("~")[-1])


# Collect directories for asex scenarios
rate_dirs = [d for d in in_dir.glob(f"SEX~0.0/GC~*") if d.is_dir()]
rate_dirs = sorted(rate_dirs, key=get_gc_rate)


for rate_dir in rate_dirs:
    gc_rate = get_gc_rate(rate_dir)
    print(f"GC rate: {gc_rate}")
    # Create a genotype directory for the given GC rate
    geno_dir = out_dir / f"SEX~0.0/GC~{gc_rate}"
    geno_dir.mkdir(parents=True, exist_ok=True)
    for ts in rate_dir.glob("*.trees"):
        geno_file = geno_dir / f"{ts.stem}.txt"
        f = open(geno_file, "w")
        ts = tskit.load(str(ts))
        # Overlay neutral mutations keeping SLiM-simulated deleterious mutations
        ts = msprime.sim_mutations(
            ts, rate=5e-9, model=msprime.SLiMMutationModel(type=1), keep=True
        )
        for var in tqdm.tqdm(ts.variants()):
            # Get the selection coefficient
            s = var.site.mutations[0].metadata["mutation_list"][0]["selection_coeff"]

            # Calculate homozygosity
            genotypes = var.genotypes
            if set(genotypes) == {0, 1} or set(genotypes) == {1}:
                genotypes = genotypes.reshape((-1, 2))
                g_scores = genotypes.sum(axis=1)
                carrier_count = np.sum(g_scores > 0)
                # Found in > 3 'libraries'
                if carrier_count < 4:
                    continue
                alt_hom_count = np.sum(g_scores == 2)
                alt_hom_prop = alt_hom_count / carrier_count

                f.write(f"{s}\t{alt_hom_prop}\n")
        f.close()
