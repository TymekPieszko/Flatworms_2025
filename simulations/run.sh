# module load Mamba/4.14.0-0

# Create envs
# snakemake --use-conda --conda-frontend mamba --conda-create-envs-only


# Dry run
# snakemake -n --cores 1 --rerun-triggers mtime

snakemake \
--executor slurm \
-j 50 \
--resources threads=16 mem_mb=160000 \
--default-resources slurm_account=biol-bdelloids slurm_partition=devel \
--group-components simulate=16 \
--use-conda --conda-frontend mamba \
--rerun-triggers mtime \
--rerun-incomplete \
--latency-wait 60000