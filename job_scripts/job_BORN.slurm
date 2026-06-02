#!/bin/bash
#SBATCH -J 600/po_STRU_md/HZO-1-3_collect_stru
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=333-08:00:00
#SBATCH --time 99-00:00:00


### Set basic var   ### MARK_slurm2pbs

JOBID="$SLURM_JOB_ID"                                  ### slurm2pbs
# 优先用 SLURM_NTASKS（更标准），否则退回 SLURM_NPROCS
NP="${SLURM_NTASKS:-${SLURM_NPROCS:-32}}"               ### slurm2pbs
NNODE=$(srun hostname | sort | uniq | wc -l)            ### slurm2pbs

# -----------------------------
# 关键修复：WORK_DIR / JOB_NAME 兜底
# -----------------------------
WORK_DIR="${WORK_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
JOB_NAME="${JOB_NAME:-${SLURM_JOB_NAME:-job}}"
mkdir -p "$WORK_DIR"

LOG_FILE="$WORK_DIR/job_${JOB_NAME}_${JOBID}.log"
HOST_FILE="$WORK_DIR/job_${JOB_NAME}_${JOBID}_${NP}c_${NNODE}n.ma"

# 预先测试日志是否可写，避免 tee 一直报错
if ! touch "$LOG_FILE" 2>/dev/null; then
  echo "[WARN] LOG_FILE 无法写入：$LOG_FILE，将只输出到标准输出，不再 tee 到文件。"
  LOG_FILE=""
fi

srun hostname | sort | uniq -c | awk '{print $2":"$1}' > "$HOST_FILE"  ### slurm2pbs

### Write basic job infomations
if [ -n "${LOG_FILE:-}" ]; then
  echo -e "The start time is: $(date +"%Y-%m-%d %H:%M:%S") \n" | tee -a "$LOG_FILE"
  echo -e "My job ID is: $JOBID \n" | tee -a "$LOG_FILE"
  echo -e "The total cores is: $NP \n" | tee -a "$LOG_FILE"
  echo -e "The hosts is: \n" | tee -a "$LOG_FILE"
  cat "$HOST_FILE" | tee -a "$LOG_FILE"
  echo -e "\n" | tee -a "$LOG_FILE"
else
  echo -e "The start time is: $(date +"%Y-%m-%d %H:%M:%S") \n"
  echo -e "My job ID is: $JOBID \n"
  echo -e "The total cores is: $NP \n"
  echo -e "The hosts is: \n"
  cat "$HOST_FILE"
  echo -e "\n"
fi

### Run APP

  # MARK_CMD  #Don't delete this line!!!
#!/bin/bash
## job script created by Gridview Jobmanager.
source ~/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
source ~/apprepo/pyatb/v1.1.2-mpi4py/scripts/env.sh
touch execution_beginning.txt
echo "begin to run task"

# =========================
# 安静 UCX / libibverbs 警告（不影响真正错误）
# =========================
export UCX_WARN_UNUSED_ENV_VARS=n
export UCX_LOG_LEVEL=error

# =========================
# 小工具：打印 +（可选）写 LOG_FILE
# =========================
log_info(){ if [ -n "${LOG_FILE:-}" ]; then echo -e "[INFO] $*" | tee -a "$LOG_FILE"; else echo -e "[INFO] $*"; fi; }
log_warn(){ if [ -n "${LOG_FILE:-}" ]; then echo -e "[WARN] $*" | tee -a "$LOG_FILE"; else echo -e "[WARN] $*"; fi; }
log_skip(){ if [ -n "${LOG_FILE:-}" ]; then echo -e "[SKIP] $*" | tee -a "$LOG_FILE"; else echo -e "[SKIP] $*"; fi; }

# 避免通配符无匹配时变成字面字符串
shopt -s nullglob

# 用 mpirun 执行，但过滤掉 UCX WARN / libibverbs Warning 行（保留其它 stderr）
mpirun_quiet() {
  # 注意：这里过滤的是“噪声警告”，其它真正的报错会保留
  mpirun "$@" 2> >(grep -v -E '(^libibverbs: Warning:|UCX  WARN )' >&2)
}

pporb()
{
case "$PP_ORB" in
  "dojo" | "DOJO" | "Dojo" | "DoJo" | "dj" )
    PP_prefix=$ABACUS_PP_DOJO
    ORB_prefix=$ABACUS_ORB_DOJO
    cp STRU STRU_bak
    python $ABACUS_DEALPP $input_stru $PP_prefix $ORB_prefix
    ;;
  "dojo-8au" )
    PP_prefix=$ABACUS_PP_DOJO
    ORB_prefix=$ABACUS_ORB_DOJO_8au
    cp STRU STRU_bak
    python $ABACUS_DEALPP $input_stru $PP_prefix $ORB_prefix
    ;;
  "dojo-auto" | "dojoauto" )
    PP_prefix=$ABACUS_PP_DOJO
    ORB_prefix=$ABACUS_ORB_DOJO_auto
    cp STRU STRU_bak
    python $ABACUS_DEALPP $input_stru $PP_prefix $ORB_prefix
    ;;
  "dojo-9au" )
    PP_prefix=$ABACUS_PP_DOJO
    ORB_prefix=$ABACUS_ORB_DOJO_9au
    cp STRU STRU_bak
    python $ABACUS_DEALPP $input_stru $PP_prefix $ORB_prefix
    ;;
  "sg15" | "sg" | "SG15"  | "SG" | "Sg15"  )
    PP_prefix=$ABACUS_PP_SG15
    ORB_prefix=$ABACUS_ORB_SG15
    cp STRU STRU_bak
    python $ABACUS_DEALPP $input_stru $PP_prefix $ORB_prefix
    ;;
  "sg15fr"  )
    PP_prefix=$ABACUS_PP_SG15_FR
    ORB_prefix=$ABACUS_ORB_SG15
    cp STRU STRU_bak
    python $ABACUS_DEALPP $input_stru $PP_prefix $ORB_prefix
    ;;
  *)
    echo "No PP_ORB set, use default STRU."
    ;;
esac
}

# ==========================================================
# 断点续算逻辑函数
# ==========================================================

# (1) 判断 SCF(abacus) 是否完成：只检查本目录一级子目录 OUT.*
#   - *.log 含 FINAL_ETOT_IS
#   - *.cube 非空
#   - 3 个 data-*.csr 都非空
is_scf_done() {
  local workdir="$1"
  local outdirs=("$workdir"/OUT.*)
  [ ${#outdirs[@]} -gt 0 ] || return 1

  local od
  for od in "${outdirs[@]}"; do
    [ -d "$od" ] || continue

    local logs=("$od"/*.log)
    [ ${#logs[@]} -gt 0 ] || continue
    local has_final=1
    local lf
    for lf in "${logs[@]}"; do
      [ -s "$lf" ] || continue
      if grep -q "FINAL_ETOT_IS" "$lf"; then
        has_final=0
        break
      fi
    done
    [ $has_final -eq 0 ] || continue

    local cubes=("$od"/*.cube)
    local has_cube=1
    local cf
    for cf in "${cubes[@]}"; do
      if [ -s "$cf" ]; then
        has_cube=0
        break
      fi
    done
    [ $has_cube -eq 0 ] || continue

    local csrs=("$od"/data-*.csr)
    [ ${#csrs[@]} -eq 3 ] || continue
    local ok_csr=0
    local sf
    for sf in "${csrs[@]}"; do
      if [ ! -s "$sf" ]; then
        ok_csr=1
        break
      fi
    done
    [ $ok_csr -eq 0 ] || continue

    return 0
  done
  return 1
}

# (4) 判断 pyatb 是否完成：
#   - Out/Polarization/polarization.dat 非空
#   - Out/running.log 含 “The calculated polarization direction is in c, P =”
is_pyatb_done() {
  local pyatb_dir="$1"
  [ -s "$pyatb_dir/Out/Polarization/polarization.dat" ] || return 1
  [ -f "$pyatb_dir/Out/running.log" ] || return 1
  grep -q "The calculated polarization direction is in c, P =" "$pyatb_dir/Out/running.log" || return 1
  return 0
}

# (2) 复制电荷密度 *.cube 到目标目录，但若目标已有 cube 则不覆盖
copy_density_to_targets() {
  local src_out_polar="$1"    # OUT.POLAR
  local parent_dir="$2"

  if [ ! -d "$src_out_polar" ]; then
    log_warn "未找到 $src_out_polar，无法复制电荷密度。"
    return 1
  fi

  local src_cubes=("$src_out_polar"/*.cube)
  if [ ${#src_cubes[@]} -eq 0 ]; then
    log_warn "$src_out_polar 下没有 *.cube，跳过复制。"
    return 0
  fi

  local d2 sub target_dir dst_dir
  for d2 in "$parent_dir"/*.*; do
    [ -d "$d2" ] || continue
    for sub in x+ y+ z+; do
      target_dir="$d2/$sub"
      if [ -d "$target_dir" ]; then
        dst_dir="$target_dir/$src_out_polar"
        mkdir -p "$dst_dir"

        local existing=("$dst_dir"/*.cube)
        if [ ${#existing[@]} -gt 0 ]; then
          log_skip "目标 $dst_dir 已存在 .cube，跳过复制以避免覆盖。"
          continue
        fi

        log_info "复制电荷密度文件到 $dst_dir"
        cp "$src_out_polar"/*.cube "$dst_dir/"
      else
        log_skip "目标子文件夹 $sub 不存在于 $d2 中，跳过。"
      fi
    done
  done
  return 0
}

# (pyatb_input + pyatb) 绑定段：
#   若 pyatb 已完成 => 两者都跳过
#   否则 => 先 pyatb_input 再 mpirun pyatb
run_pyatb_segment() {
  local workdir="$1"   # 当前计算目录（含 INPUT-scf 的那个目录）

  if [ -d "$workdir/pyatb" ] && is_pyatb_done "$workdir/pyatb"; then
    log_skip "$workdir/pyatb: 检测到 polarization.dat + running.log 字段齐全，跳过 pyatb_input + pyatb。"
    return 0
  fi

  log_info "$workdir: pyatb 未完成，执行 pyatb_input --polar --mp=0.08"
  ( cd "$workdir" && pyatb_input --polar --mp=0.08 )

  if [ -d "$workdir/pyatb" ]; then
    log_info "$workdir/pyatb: 开始运行 pyatb（mpirun）"
    ( cd "$workdir/pyatb" && mpirun_quiet -genv OMP_NUM_THREADS=1 -np "$NP" pyatb )
  else
    log_warn "$workdir: pyatb_input 执行后仍未生成 pyatb 目录，跳过 pyatb。"
  fi
}

# ==========================================================
# 主流程
# ==========================================================

ROOT_DIR="$PWD"
export OMP_NUM_THREADS=1

# ---------------------------
# 先在 0.no-move 目录执行
# ---------------------------
if [ -d "$ROOT_DIR/0.no-move" ]; then
  cd "$ROOT_DIR/0.no-move" || exit 1

  if [ ! -f INPUT-scf ]; then
    log_warn "不在POLAR计算文件夹内(缺少 INPUT-scf)，退出"
    exit 1
  fi

  # (1) SCF 断点续算判断
  if is_scf_done "$PWD"; then
    log_skip "0.no-move: 检测到 SCF 已完成，跳过 abacus。"
  else
    log_info "0.no-move: SCF 未完成，执行 abacus（mpirun）..."
    cp INPUT-scf INPUT
    mpirun_quiet -np "$NP" abacus > scf.out
  fi

  log_info "当前目录是 '0.no-move'：复制电荷密度并进行 pyatb 相关流程。"

  out_polar_dir="OUT.POLAR"
  parent_dir="$(dirname "$PWD")"

  # (2) 无论 scf 是跑完还是跳过，都尝试复制，但不覆盖已有 cube
  copy_density_to_targets "$out_polar_dir" "$parent_dir" || log_warn "电荷密度复制步骤有问题（已提示），继续。"

  # (3&4) pyatb_input + pyatb 绑定段
  run_pyatb_segment "$PWD"
fi

# ---------------------------
# 再遍历 ROOT_DIR 下其他目录（x+/y+/z+）
# ---------------------------
cd "$ROOT_DIR" || exit 1

for dir in "$ROOT_DIR"/*; do
  if [ -d "$dir" ] && [ "$(basename "$dir")" != "lost+found" ] && [ "$(basename "$dir")" != "0.no-move" ]; then
    for sub in x+ y+ z+; do
      target_dir="$dir/$sub"
      if [ -d "$target_dir" ]; then
        cd "$target_dir" || continue

        if [ ! -f INPUT-scf ]; then
          log_warn "$target_dir: 缺少 INPUT-scf，退出"
          exit 1
        fi

        # (1) SCF 断点续算判断
        if is_scf_done "$PWD"; then
          log_skip "$target_dir: 检测到 SCF 已完成，跳过 abacus。"
        else
          log_info "$target_dir: SCF 未完成，执行 abacus（mpirun）..."
          cp INPUT-scf INPUT
          mpirun_quiet -np "$NP" abacus > scf.out
        fi

        # (3&4) pyatb_input + pyatb 绑定段
        run_pyatb_segment "$PWD"

        cd "$ROOT_DIR" || exit 1
      fi
    done
  fi
done

cd "$ROOT_DIR" || exit 1
zstar deal
echo "End Time: $(date)"

  # MARK_BASH #Don't delete this line!!!

if [ $? -ne 0 ]; then
  exit 1
fi

if [ -n "${LOG_FILE:-}" ]; then
  echo "The end time is: $(date +"%Y-%m-%d %H:%M:%S")" | tee -a "$LOG_FILE"
else
  echo "The end time is: $(date +"%Y-%m-%d %H:%M:%S")"
fi
