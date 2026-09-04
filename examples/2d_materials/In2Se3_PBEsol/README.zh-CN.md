# Alpha-In2Se3：PBEsol BEC 文献对照

本目录对应论文逐层 BEC 表中的铁电单层 PBEsol 数据，不是旁边
`In2Se3` 目录的旧 PBE+D3(0) 示例，也不是 `Shared_Response` 中新优化的
PBE 效率基准。

实际输入采用 PBEsol、D3(0)、100 Ry、SCF 阈值 1e-8 和
`kspacing 0.1 0.1 1`；附带 In/Se ONCV 赝势和 10 au 轨道。
归档中的显式 D3 参数保持原样。

安装 ZStar、PYATB 并配置 ABACUS 路径及 MPI/OMP 后，执行 `bash run.sh`。
脚本将 `run/` 的干净输入复制到 `work/`，按原始笛卡尔单边方案生成
15 个位移和一个参考结构，串行完成计算并后处理；原始结果不被覆盖。

`results/` 保留输入、结构、SCF 日志、极化输出和原始/对称化张量。
`response_observations.json` 保留实际位移、偶极变化和源文件哈希，可离线
重建张量。大型矩阵和电荷 cube 不随案例分发，可通过重新运行 SCF 生成。
该案例只验证 BEC，不声明此结构的声子稳定性或静态介电响应。

沿 +z 的层序为 Se(1)、In(1)、Se(c)、In(2)、Se(2)。旧版带原子编号的
ZStar 张量表是位移优先排列；Phonopy 的 BORN 是极化优先排列。
