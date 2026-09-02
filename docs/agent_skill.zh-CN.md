# 面向智能体的 ZStar 工作流

ZStar 随软件发行符合规范的 Agent Skill：`run-zstar-workflows`。它将命令行接口、
物理约定、断点续算状态和输出完成判据组织为可复用指令，使兼容的编程或科研智能体
能够自动发现并正确调用 ZStar。

该 Skill 不替代科研判断，也不会自动获得启动计算的权限。它帮助智能体选择正确的
工作流、检查运行条件、保留来源信息、复用已完成阶段，并依据可观察产物判断软件任务
是否完成。

## 命名与目录规范

Skill 采用可移植的 Agent Skills 目录结构：

```text
run-zstar-workflows/
|-- SKILL.md
|-- agents/openai.yaml
|-- references/
`-- scripts/preflight.py
```

目录名与 YAML 中的 `name` 完全一致，均为 `run-zstar-workflows`：小写 ASCII、
连字符分隔、动作导向，并且少于 64 个字符。

## 安装

安装 ZStar 后，可将随 wheel 提供的 Skill 安装到默认 Codex 技能目录：

```bash
pip install -U zstar
zstar skill install
```

安装后新建一个智能体会话，使其重新发现技能。升级 ZStar 后可覆盖旧 Skill：

```bash
zstar skill install --force
```

也可以安装到其他兼容智能体框架的技能父目录：

```bash
zstar skill install --dest /path/to/skills
zstar skill path
```

源码仓库中的 Skill 位于 `zstar/agent_skills/run-zstar-workflows/`，兼容的 Skill
安装器也可以直接从该仓库路径安装。

## 调用

显式调用采用标准技能名，例如：

```text
使用 $run-zstar-workflows 检查这个 BaTiO3 bulk BEC 工作区，执行 preflight，
生成 Slurm dry-run 脚本，但不要提交任务。
```

当请求明确涉及 ZStar 计算的准备、执行或验证时，技能描述也支持自动发现。

## 机器可读的 preflight

智能体在构造或启动工作流之前，应先检查工作区：

```bash
zstar skill preflight --root . --lane bec --dim bulk
zstar skill preflight --root . --lane raman --dim 2d
zstar skill preflight --root . --lane ir --dim molecule
zstar skill preflight --root . --lane database --dim 1d
```

该命令只向标准输出写 JSON，不修改工作区。报告包括：

- ZStar 和 Python 版本；
- 科学任务类型与维度约定；
- 必需输入阻塞项和环境警告；
- 外部程序与已有产物；
- `.zstar/stages/*.json` 的状态统计和失败诊断；
- 根据必需输入给出的 `ready` 布尔值，但不会伪装成物理收敛认证。

## 智能体契约

Skill 固化了不应随语言模型变化的关键约束：

| 契约 | 智能体行为 |
| --- | --- |
| 参考态优先 | 位移 BEC 之前完成 `0.no-move` 并通过绝缘性门控。 |
| 低维分流 | 一维纳米线与二维薄膜沿周期轴使用 Berry 响应，沿开放轴使用 cube 偶极积分。 |
| 一维边界 | 可执行沿 `z` 周期的 ABACUS + PYATB BEC 与 Gamma 点光谱；没有真正 1D Coulomb cutoff 时，不得宣称有限波矢极性声子已经完成。 |
| 断点续算 | 复用 `.zstar` 状态并重复同一串行执行命令。 |
| Raman 差分 | 使用正、负简正坐标位移。 |
| 完成判据 | 检查指定产物、状态记录、单位和物理约定。 |
| 权限边界 | 不把准备脚本理解为允许提交或启动昂贵远程计算。 |

因此，这个接口由文本技能、确定性 JSON 检查、稳定 CLI 和机器可验证产物共同组成，
不依赖某一个特定语言模型也可以阅读和使用。
