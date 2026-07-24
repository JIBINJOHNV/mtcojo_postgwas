"""
postgwas_runner.py - PostGWAS Docker Harmonisation Launcher for mtcojo_postgwas

Merges GCTA mtCOJO results with in-memory VCF coordinates (0 VCF re-reads),
generates /mnt/disks/sdd/GWAS_sumstat2/gwas2vcf_input2.tsv manifest, and runs Docker harmonisation.
"""

import os
import shutil
import subprocess
import polars as pl

def run_postgwas_harmonisation(
    cma_file: str,
    target_coords: pl.DataFrame,
    out_prefix: str,
    sdd_dir: str,
    defaults_yaml: str,
    resource_folder: str,
    docker_image: str = "jibinjv/postgwas:1.4",
    docker_platform: str = "linux/amd64",
    nthreads: int = 23,
    max_mem: str = "50G",
    liftover: str = "No"
) -> None:
    """
    Execute PostGWAS Docker harmonisation.

    Args:
        cma_file (str): Path to GCTA mtCOJO .mtcojo.cma output file.
        target_coords (pl.DataFrame): In-memory DataFrame containing [SNP, CHR, POS, SI].
        out_prefix (str): mtCOJO output prefix name.
        sdd_dir (str): Host directory mounted as /mnt/disks/sdd/.
        defaults_yaml (str): Path to harmonisation.yaml config file.
        resource_folder (str): Path to gwas2vcf reference resources.
        docker_image (str): Docker image tag.
        docker_platform (str): Docker platform flag.
        nthreads (int): Thread count for Docker container.
        max_mem (str): Max memory allocation for Docker container.
        liftover (str): "Yes" or "No".
    """
    print("\n" + "=" * 70)
    print(" [PostGWAS Runner] Starting PostGWAS Docker Harmonisation")
    print("=" * 70)

    sdd_dir = os.path.abspath(sdd_dir)
    gwas_dir = os.path.join(sdd_dir, "GWAS_sumstat2")
    out_folder = os.path.join(sdd_dir, "03_harmonised_output")
    os.makedirs(gwas_dir, exist_ok=True)
    os.makedirs(out_folder, exist_ok=True)

    gwas_name = os.path.basename(out_prefix)
    print(f"[PostGWAS Runner] Merging mtCOJO results with in-memory coordinates (0 VCF re-reads)...")
    
    df_cma = pl.read_csv(cma_file, separator="\t", truncate_ragged_lines=True)
    df_merged = df_cma.join(target_coords, on="SNP", how="left")

    sumstat_tsv = os.path.join(gwas_dir, f"{gwas_name}_harmonised_sumstat.tsv")
    df_merged.select([
        pl.col("CHR").fill_null("1").alias("CHR"),
        pl.col("POS").fill_null(0).alias("POS"),
        pl.col("SNP"), pl.col("A1"), pl.col("A2"), pl.col("freq"),
        pl.col("bC"), pl.col("bC_se"), pl.col("bC_pval"),
        pl.col("N"), pl.col("SI")
    ]).write_csv(sumstat_tsv, separator="\t")
    print(f"[PostGWAS Runner] Enriched summary stats written to: {sumstat_tsv}")

    # Generate 25-column config manifest gwas2vcf_input2.tsv
    config_path = os.path.join(gwas_dir, "gwas2vcf_input2.tsv")
    resource_parent = os.path.dirname(os.path.abspath(resource_folder.rstrip("/")))
    resource_basename = os.path.basename(os.path.abspath(resource_folder.rstrip("/")))
    container_resource_root = "/opt/postgwas_resources"
    container_resource_folder = f"{container_resource_root}/{resource_basename}"

    manifest_data = [{
        "sumstat_file": f"/mnt/disks/sdd/GWAS_sumstat2/{gwas_name}_harmonised_sumstat.tsv",
        "gwas_outputname": f"{gwas_name}_harmonised",
        "chr_col": "CHR", "pos_col": "POS", "snp_id_col": "SNP",
        "ea_col": "A1", "oa_col": "A2", "eaf_col": "freq",
        "beta_or_col": "bC", "se_col": "bC_se", "imp_z_col": "NA", "pval_col": "bC_pval",
        "ncontrol_col": "N", "ncase_col": "NA", "ncontrol": "NA", "ncase": "NA",
        "imp_info_col": "SI", "infofile": "NA", "infocolumn": "NA", "eaffile": "NA", "eafcolumn": "NA",
        "liftover": liftover, "chr_pos_col": "NA",
        "resource_folder": container_resource_folder, "resourse_folder": container_resource_folder,
        "output_folder": "/mnt/disks/sdd/03_harmonised_output"
    }]
    pl.DataFrame(manifest_data).write_csv(config_path, separator="\t")
    print(f"[PostGWAS Runner] Generated config manifest: {config_path}")

    # Copy defaults YAML config
    dest_yaml_path = os.path.join(sdd_dir, "harmonisation.yaml")
    if os.path.exists(defaults_yaml):
        shutil.copyfile(defaults_yaml, dest_yaml_path)
        print(f"[PostGWAS Runner] Copied defaults config from: {defaults_yaml}")
    else:
        raise FileNotFoundError(f"Defaults YAML config file not found: {defaults_yaml}")

    docker_cmd = [
        "docker", "run", f"--platform={docker_platform}",
        "-u", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{sdd_dir}:/mnt/disks/sdd/",
        "-v", f"{resource_parent}:{container_resource_root}:ro",
        docker_image, "postgwas", "harmonisation",
        "--nthreads", str(nthreads), "--max-mem", max_mem,
        "--config", "/mnt/disks/sdd/GWAS_sumstat2/gwas2vcf_input2.tsv",
        "--defaults", "/mnt/disks/sdd/harmonisation.yaml"
    ]
    
    print(f"[PostGWAS Runner] Executing Docker Command:\n{' '.join(docker_cmd)}\n")
    subprocess.run(docker_cmd, check=True)
    print("[PostGWAS Runner] PostGWAS Harmonisation Completed Successfully!")
