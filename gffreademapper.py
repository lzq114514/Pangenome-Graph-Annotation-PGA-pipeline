#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sys

# -------------------
# 参数检查
# -------------------
if len(sys.argv) != 2:
    print("❌ Error: Please provide the root directory path.")
    print("Usage: python script.py <input_dir>")
    sys.exit(1)

input_dir = sys.argv[1]
if not os.path.exists(input_dir):
    print(f"❌ Error: Input directory {input_dir} does not exist!")
    sys.exit(1)

# workflow1 文件夹
output_dir = os.path.join(input_dir + 'workflow1')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# -------------------
# 遍历子目录，生成 .pep
# -------------------
def get_all_subdirectories(root_dir):
    subdirs = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if dirpath.startswith(output_dir):
            continue
        subdirs.append(dirpath)
    return subdirs

all_subdirs = get_all_subdirectories(input_dir)

for sub_dir_path in all_subdirs:
    fna_file, gff_file = None, None
    for file in os.listdir(sub_dir_path):
        if file.endswith(".fna"):
            fna_file = os.path.join(sub_dir_path, file)
        elif file.endswith(".gff"):
            gff_file = os.path.join(sub_dir_path, file)

    if fna_file and gff_file:
        gff_prefix = os.path.splitext(os.path.basename(gff_file))[0]
        target_file = os.path.join(output_dir, f"{gff_prefix}.pep")
        command = f"gffread {gff_file} -g {fna_file} -y {target_file}"
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"✅ gffread OK: {target_file}")
        except subprocess.CalledProcessError as e:
            print(f"❌ gffread failed: {sub_dir_path}: {e}")
    else:
        if os.listdir(sub_dir_path):
            print(f"⚠️ Missing fna/gff in {sub_dir_path}")

# -------------------
# 清理 .pep 文件，直接覆盖
# -------------------
for file in os.listdir(output_dir):
    file_path = os.path.join(output_dir, file)
    if os.path.isfile(file_path) and file.endswith(".pep"):
        print(f"🔧 Cleaning {file}...")
        try:
            with open(file_path, 'r') as f_in:
                lines = f_in.readlines()
            with open(file_path, 'w') as f_out:
                for line in lines:
                    if line.startswith('>'):
                        f_out.write(line)
                    else:
                        f_out.write(line.replace('.', ''))
            print(f"  ✅ Cleaned {file} (dots removed)")
        except Exception as e:
            print(f"❌ Error cleaning {file}: {e}")

# -------------------
# 生成 SLURM 脚本提交 emapper
# -------------------
for file in os.listdir(output_dir):
    file_path = os.path.join(output_dir, file)
    if os.path.isfile(file_path) and file.endswith(".pep"):
        emapper_output_file = os.path.join(output_dir, f"{os.path.splitext(file)[0]}_emapper")
        slurm_script = os.path.join(output_dir, f"run_emapper_{os.path.splitext(file)[0]}.sh")
        try:
            with open(slurm_script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"#SBATCH --job-name=emapper_{file}\n")
                f.write("#SBATCH --partition=hebhcnormal01\n")
                f.write("#SBATCH --nodes=1\n")
                f.write("#SBATCH --ntasks-per-node=20\n")
                f.write("#SBATCH --error=%j.err\n")
                f.write("#SBATCH --output=%j.out\n\n")
                f.write(f"emapper.py -i {file_path} -o {emapper_output_file} "
                        f"--cpu 10 --itype proteins -m diamond --override\n")
            subprocess.run(f"sbatch {slurm_script}", shell=True, check=True)
            print(f"  ✅ SLURM job submitted for {file}")
        except subprocess.CalledProcessError as e:
            print(f"❌ sbatch submission failed: {file}: {e}")
        except Exception as e:
            print(f"❌ Unexpected error: {file}: {e}")
