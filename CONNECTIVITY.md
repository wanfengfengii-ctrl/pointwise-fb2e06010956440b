# 周期像连通成果（新增公开契约）

现有正常原胞闭环、输入数值及 V1 几何成果见 README.md。本次只新增周期像连通加工；
不生成超胞、不改变晶格、不猜键，不要求最短路径、展开走向或枚举全部证据。
生产实现仅用 Python 标准库，继续使用已有 Archive、Tasks 和 Results。

## 调用及有限边界

`Tasks(archive).run_connectivity(product_id, revision, queries)` 返回持久成果整数 ID。
产品、修订及成果 ID 的合法性和缺失语义沿用 V1。所引用修订至多 16 个原子、32 条键；
仅连通加工检查此上限，不限制原来的接纳、几何加工或历史成果。
queries 必须为长度 0～16 的 list，每项是恰含 u、v、translation 的 dict。
u、v 为该修订中存在的原子 ID；translation 是三个严格 Python int 构成的 list 或 tuple，
不接受 bool、Fraction、浮点或有理数对。整数大小不另作限制。输出 query 的 translation 统一为 list。
空查询合法；空原胞仅可接收空查询。重复查询逐项保留，不合并、不排序。
不增加输入文件格式；已有 cell 格式、单位、原子及键身份完全沿用 V1。

原键 (id,u,v,t) 在每个 n∈Z³ 处连接 (u,n) 与 (v,n+t)。
本分析允许正反向行走；反向仍用同一原键 ID，贡献相反偏移和负使用次数，
不改 V1 有向位移、长度平方或键身份语义。
逐项准确判断 (u,0) 到 (v,translation) 是否存在有限走向。
允许全量重算与精确整数消元；不能以实数/有理数张成、浮点容差或固定展开窗口代替整数可达性。

## 读取和整数链证据

`Results(archive).read(result_id)` 对新类型返回恰含 source、answers 的独立字典。
source 与 V1 一样恰含 product_id、revision、digest、calculation，
calculation 为 `connectivity-v1`，digest 仍是完整规范 cell 的既有摘要。
answers 是按输入顺序排列的 list；每项恰含 query、reachable、coefficients：
query 为独立的规范查询字典，reachable 为 bool。
不可达时 coefficients 为 None；可达时为原键 ID 到非零严格整数系数的 dict，
不记零项，没有净使用的键时用空字典。证据任选，不要求最短、唯一或最小系数。

系数是原键正向净使用次数，必须满足以下全部条件：
- u、v 在原胞有限无向图的同一连通分量，孤立原子自身也是分量；
- 每个出现的键 ID 都属于该分量，不能借用其他分量的闭合连接；
- 对每个原子，系数的流入减流出等于该原子是 v 的指示量减去它是 u 的指示量；
  自像键在端点收支上相消；
- 所有键 offset 乘系数后的三维整数和恰为 translation。

这是整数链证据，不要求逐步走向。仅有一个原子 u 和键 a=(u,u,(2,0,0)) 时，
平移 (1,0,0) 不可达，(2,0,0) 可达；增加独立键 b=(u,u,(3,0,0)) 后，
平移 (1,0,0) 可达，证据可为 {"a":-1,"b":1}。
若 u、v 之间仅有偏移 (1,0,0)、(3,0,0) 两条键，则 u 到 v 只能到达奇数 x 平移，
(2,0,0) 不可达，即使对应实数方程有解。

## 固定成果和离库交付

源修订、摘要、完整查询及结论和证据整体持久保存，Results.read 不得重新求解。
更正来源、修改原查询或返回容器、关闭重开均不得改变旧成果。
新类型沿同一 Results.export 返回 JSON 文本；包恰含 format、version、source、result，
format 为 `crystal-connectivity-evidence`、version 为整数 1。
source 沿用 V1 导出中的 product_id、revision、digest、完整规范 cell 和 raw_cell；
result 与新类型 read 返回的完整字典相同，布尔、整数和 None 用原生 JSON 表示。
导出保留原始与规范来源的一致性、摘要及头体绑定核验；不调用连通求解器重新生成结论。
离库核账只需这个包，不打开数据库。V1 导出格式及既有值相等语义不变。

非法查询、端点或本次规模超限抛 ValueError，整次不新增任务或成果；
合法但缺失的产品、修订、成果抛 KeyError。
计算或发布异常原样传播，不留本次任务或部分成果，旧成果不变。
保留全部 V1 API 和原测试；已建 V1 库及其旧修订、几何成果可继续读取并接纳新任务。
可扩展现有存储或增加同库结果表示，不指定内部表、值对象或求解组织。
不增加并发调度、增量索引或崩溃恢复协议。

## 样例、运行及交付

connectivity_cases.py 是非答案材料：只有输入、手算结论和独立正例证据核账，
没有整数求解器或目标生产实现。可先运行
`python3.11 -B connectivity_cases.py --self-check`，只核验夹具输入和证据检查器。原有 verify.py、demo.py 及运行清单仍可执行 V1。

新增能力完成后：
`python3.11 -B connectivity_cases.py --memory` 检查普通连通、非平凡整数格、
身份/方向、输入拒绝、发布拒绝和固定读取；
`python3.11 -B connectivity_cases.py` 再做真实磁盘关闭重开及离库核账，并展示
偏移 2 的旧修订拒绝平移 1、加入偏移 3 的新修订接纳平移 1。
断言或业务异常使该命令以非零码退出。脚本不能替代交付者的独立核账与正常回归。

保留全部原测试字节，提供简短的表示、原键证据及固定成果交接不变量说明；
补充独立小例核账与必要的失败、重开验证。
更新现有 .pointwise-runtime.json，让 verify 核验完整 V1 和新链，demo 展示上述修订、
持久读取及离库证据。保留现有清单约定：version=1、kind=cli，setup 是 argv 数组列表，
verify、demo 各是可运行的 argv 数组；无需安装依赖时 setup 可仍为空。
