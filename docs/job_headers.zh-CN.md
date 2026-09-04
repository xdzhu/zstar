# 作业 header 与计算环境

[English tutorial](job_headers.md)

ZStar 为计算目录生成一个串行、可断点续算的脚本。调度队列、资源和
`module`/环境加载命令放在 header；软件路径、MPI/OMP、赝势和轨道目录放在配置中。
生成脚本不会自动提交任务。

## 三个位置，按优先级选一个

| 优先级 | 名称 | 位置 |
| --- | --- | --- |
| 1 | Specified | `--header /path/to/header.sh` |
| 2 | Current | 工作流根目录中的 `header.sh`，即 `--root`（默认当前目录） |
| 3 | Global | `~/.zstar/header.sh` |
| 4 | Default | 都不存在时，生成可直接编辑的默认模板 |

不合并 header，不向上搜索其他目录。特别指定的文件不存在时直接报错；
选中的文件为空也报错。确实不需要任何内容时，可使用只有一行注释的文件。
选中 header 后，以它的调度资源为准，不再混入默认的队列、节点、账号和时限。

先设置计算器及并行参数：

```bash
zstar config set executables.abacus /opt/abacus/bin/abacus --user
zstar config set execution.mpi 1 --user
zstar config set execution.omp 40 --user
zstar config set abacus.pseudo_dir /data/PSEUDO --user
zstar config set abacus.orbital_dir /data/ORBITAL --user
```

去掉 `--user` 则写入工作流本地配置。旧配置 `execution.tasks` 和
`execution.cpus_per_task` 继续可用；同时存在时，以 `mpi` 和 `omp` 为准。
显式的 `--tasks`、`--cpus-per-task` 可覆盖配置。
并行设置须与 header 申请的资源匹配；软件不从任意 shell 命令反推 MPI/OMP。

## Slurm 示例

在计算根目录的 `header.sh` 中填写以下内容，或者放到全局位置
`~/.zstar/header.sh`。队列、模块和路径须替换成实际集群配置。

```bash
#!/usr/bin/env bash
#SBATCH --job-name=zstar
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00
#SBATCH --output=zstar-%j.out

module load compiler mpi
source /opt/conda/etc/profile.d/conda.sh
conda activate zstar
```

调度指令必须写在所有 shell 命令前面。然后在计算目录执行：

```bash
zstar bec pre --stru STRU
zstar bec job --system slurm
sbatch run_zstar_born.slurm
```

切换集群时可特别指定：

```bash
zstar bec job --system slurm --header /data/cluster-b/header.sh
```

`zstar phonon job`、`zstar spectra job` 及各计算器的脚本生成入口使用相同规则。

## Torque 与本地运行

Torque 的同等资源示例：

```bash
#!/usr/bin/env bash
#PBS -N zstar
#PBS -q compute
#PBS -l nodes=1:ppn=40
#PBS -l walltime=24:00:00

module load compiler mpi
source /opt/conda/etc/profile.d/conda.sh
conda activate zstar
```

```bash
zstar bec job --system torque
qsub run_zstar_born.pbs
```

本地 `--system shell` 的 header 只保留环境设置，不放 `#SBATCH` / `#PBS`。
生成后执行 `bash run_zstar_born.sh`。header 调度类型与 `--system` 不一致会报错，
避免全局 Slurm header 被误用；此时用 `--header` 指向正确文件。

## 复现与断点续算

header 内容直接嵌入生成脚本；路径、SHA-256、选中层级写入
`.zstar/job_header.json`。以后修改原 header 不会偷偷改变已生成的脚本，
需要重新生成。脚本采用 Bash 严格错误检查，从工作流根目录执行环境命令。

提交前检查资源、账号、模块及 MPI/OMP 的一致性。重新运行相同脚本会根据
`.zstar` 中的阶段记录继续未完成计算；`zstar bec stat` 可查看进度。
这些改动不改变参考结构优先计算和电荷密度独立复制的行为。
旧 `--env-script` 仍兼容，新工作流建议将环境命令直接放到 header 中。

三个位置都没有文件时，仍生成完整默认模板；顶部简短注释指出三个位置并链接
完整教程。一次性任务可直接修改生成脚本，常用环境则保存为 header。
Phonopy、PYATB 位于激活的 Python 环境中；DFT 软件使用 `PATH` 或已配置的路径。
