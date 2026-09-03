# ZStar 命令行参考

ZStar 按科学对象组织公共命令。凡是需要实际计算的工作流，都尽量采用同一套生命
周期：

```text
pre -> run -> job -> stat -> post
```

- `pre` 生成输入，并写入 `.zstar/<family>.json` 工作流清单。
- `run` 串行执行已准备的阶段，自动跳过已完成阶段并支持断点续算。
- `job` 为整条串行链生成一个 shell、Slurm 或 Torque 驱动脚本。
- `stat` 汇总准备、运行、完成、失败和阻塞状态。
- `post` 校验并收集科学结果。

`prepare`、`status`、`collect`/`deal` 和 `script` 仍作为兼容别名可用；新文档与
自动化应使用上面的规范短动词。

## 公共命令树

| 功能族 | 规范动作 | 用途 |
| --- | --- | --- |
| `zstar bec` | `pre/job/run/stat/post` | 极化、APT/BEC 与 `BORN`；支持 ABACUS + PYATB、VASP、CP2K、QE。 |
| `zstar phonon` (`ph`) | `pre/job/run/stat/post/irrep` | 位移、串行力计算、力常数、频率与 Gamma 点不可约表示。 |
| `zstar spectra` | `pre/job/run/stat/post` | ABACUS + PYATB、VASP、CP2K、QE 的 IR 与 Raman 工作流。 |
| `zstar dielectric` (`diel`) | `static` (`zero`)、`freq`、`optics` | 晶格静态响应、频率相关振动响应与电子光学响应。 |
| `zstar backend list` | `--check`、`--json`、`--discover` | 列出已实现能力，并可检查本机程序或第三方插件。 |
| `zstar config` | `init/show/set/check` | 分层管理计算软件路径与运行配置。 |
| `zstar response` | `validate/import-bec/import-abacus/import-phonopy/intrinsic` | 计算器无关响应文档及低维本征响应。 |
| `zstar density` | `vasp-cube/qe-input/qe-sidecar/cp2k-block/sidecar` | 电荷密度导出适配器与来源 sidecar。 |
| `zstar stru` | `convert/wyckoff` | 结构转换与对称性检查。 |
| `zstar data` | `db/qnep` | 可追溯 BEC/High-K 数据库和 qNEP 训练数据。 |
| `zstar skill` | `install/path/preflight` | 安装或检查随包发布的 Agent Skill，并执行只读预检查。 |
| `zstar pot` | 选项驱动 | 轴向曲线、平面图、方向曲线、真空势差和镜面对称破缺。 |

旧的细粒度命令（`gen`、`workflow`、`deal`、`postph`、`ir`、`raman` 及各类
`*-bec`）继续作为兼容层和专家接口保留。后端公共入口只有
`zstar backend list`。

## 计算软件路径配置

在项目内初始化配置并一次性设置程序路径：

```bash
zstar config init
zstar config set executables.abacus /opt/abacus/bin/abacus
zstar config set executables.pyatb /opt/pyatb/bin/pyatb
zstar config set executables.vasp /opt/vasp/bin/vasp_std
zstar config set executables.cp2k /opt/cp2k/bin/cp2k.psmp
zstar config check
```

可配置键为 `abacus`、`pyatb`、`vasp`、`cp2k`、`qe_pw`、`qe_ph`、
`qe_dynmat` 和 `phonopy`。优先级从低到高依次为：内置默认值、用户配置、项目
配置、环境变量；单次命令行覆盖的优先级最高。项目配置位于
`.zstar/config.toml`；Linux 用户配置默认位于 `~/.config/zstar/config.toml`，
Windows 位于 `%APPDATA%/zstar/config.toml`。`ZSTAR_CONFIG` 可指定其他用户
配置文件。

环境变量使用 `ZSTAR_ABACUS_EXECUTABLE`、`ZSTAR_QE_PW_EXECUTABLE` 等名称。
配置文件保存计算软件可执行文件及 ABACUS 的全局资源目录；ZStar 会为 shell/Torque 添加
`mpirun -np N`，为 Slurm 添加 `srun --ntasks=N`。模块加载、conda 激活等环境
初始化写在 `job --env-script` 指定的脚本中。

## 集群资源与运行环境

队列和资源参数属于某一次具体任务，因此在生成作业脚本时指定，不写入计算软件
路径配置：

```bash
zstar bec job --root . --system slurm \
  --nodes 1 --tasks 28 --cpus-per-task 1 \
  --queue compute --account PROJECT --walltime 24:00:00 \
  --env-script ./env.sh
```

`zstar phonon job` 和 `zstar spectra job` 也支持同样的资源选项。直接在本地节点
运行使用 `--system shell`；Torque/PBS 使用 `--system torque` 或 `--system pbs`。
生成的脚本包含相应的调度器指令，Slurm 使用 `srun`，shell/Torque 使用 `mpirun`。
`env.sh` 中放置 `module load`、`conda activate` 及集群特有的环境变量。使用
`--dry-run` 可以只检查生成脚本而不启动计算器。

## ABACUS 赝势与数值轨道

`zstar bec pre` 会在创建位移任务之前自动解析 ABACUS 的 `.upf` 和 `.orb`
文件。对于某个案例需要覆盖默认目录时，可以直接指定：

```bash
zstar bec pre --stru STRU \
  --pp /path/to/PSEUDO \
  --orb /path/to/ORBITAL
```

常用目录可以写入用户全局配置：

```bash
zstar config set abacus.pseudo_dir /data/PSEUDO --user
zstar config set abacus.orbital_dir /data/ORBITAL --user
```

ZStar 首先按照 `STRU` 中给出的精确文件名查找；精确文件名不存在时，只有在
某元素的前缀匹配结果唯一时才会自动采用，例如 `Si_*.upf` 或 `Si_*.orb`。
如果结果不唯一，命令会列出全部候选文件，并提示用户填写精确文件名或指定更
窄的目录，绝不会静默选择第一个文件。

原始 `STRU` 不会被修改。解析后的副本保存在 `.zstar/STRU.resolved`，实际选用
的文件及其 SHA256 校验值保存在 `.zstar/assets.json`，并会被复制到每个生成的
ABACUS 任务目录中。

## 代表性生命周期

以下示例均在已经放好计算器输入文件的工作目录中执行。ABACUS + PYATB 三维 BEC：

```bash
zstar bec pre --stru STRU
zstar bec job --root . --system slurm --tasks 28 --env-script env.sh
zstar bec stat --root .
zstar bec post --root .
```

`zstar bec pre` 默认使用 ABACUS + PYATB 和 forward 前向有限差分，因此规范示例
省略计算器和 `--pyatb` 开关。只有切换后端时才指定 `--calculator cp2k`、`vasp`
或 `qe`。

声子与模式分类：

```bash
zstar phonon pre --stru STRU --dim "2 2 2"
zstar phonon job --system slurm --tasks 28
zstar phonon stat
zstar phonon post
zstar phonon irrep
```

IR 与 Raman：

```bash
zstar spectra pre --stru STRU
zstar spectra job --system slurm --tasks 28
zstar spectra stat
zstar spectra post
```

静态与频率相关介电响应：

```bash
zstar dielectric static
zstar dielectric freq
```

## 静电势能力闭环

`zstar pot` 保留了 ZStar 论文中展示的完整分析能力：

```bash
zstar pot --cube ElecStaticPot.cube --axes z \
  --plane xy --plane-average --tile 5 5 \
  --direction a+b --mirror-test \
  --vacuum-sides --vacuum-window 0.75 --outdir potential
```

该命令输出一维平均势、带中心虚线单胞框的扩胞平面势图、垂直于指定方向的平面
平均曲线、单周期最佳镜面中心与非对称度，以及上下表面的真空势差。镜面指标评价
的是选定周期方向内部相对于完美镜面对称轮廓的破缺，而不是比较 `a+b` 与 `a-b`。

## 兼容命令对照

| 规范命令 | 仍可使用的旧写法 |
| --- | --- |
| `zstar bec pre` | `zstar gen` |
| `zstar bec run/job/stat` | `zstar workflow run/status/script` |
| `zstar bec post` | `zstar deal` |
| `zstar phonon pre/post/irrep` | `zstar ph/postph/irrep` |
| `zstar spectra pre/job/run/stat/post` | `prepare/script/run/status/collect`；底层 `ir`、`raman` 继续可用 |
| `zstar dielectric static/freq` | `zstar calc/freq` |
| `zstar stru convert/wyckoff` | `zstar vasp/wyckoff` |
| `zstar skill install/path/preflight` | `zstar agent-skill ...` |
