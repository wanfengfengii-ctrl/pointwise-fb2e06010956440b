# 晶格材料实验室原胞加工账本

Python 3.11，仅标准库。实验记录整理人员提供晶格、原子及已确认的周期键；
加工员登记完整送件和更正，下游模拟准备人员读取精确几何成果或离库核账。
本软件处理原胞，不猜键、不取最短像、不进行化学或结构优化。

## 一天的加工流程

1. 接收产品ID和完整送件，为原始容器建立独立快照。
2. 把晶格长度精确换算为nm，保留原单位审计材料。
3. 核验可逆晶格和半开坐标，建立原子身份索引。
4. 核验周期键端点、身份和整数偏移，提交新不可变修订及当前指针。
5. 下游按明确产品和修订固定任务，任务不追随当前指针。
6. 计算体积、原子位置和逐键账，核对身份覆盖与物理值。
7. 在同一事务提交任务、成果头和成果体，下游只读取完整成果。
8. 保存更正并分别计算新旧版本，关闭重开，导出旧版证据供离库复算。

## 运行

从项目根目录执行：

```sh
python3.11 -B verify.py
python3.11 -B demo.py
```

verify发现并执行全部tests/test_*.py（含连通链），随后执行
connectivity_cases.py --memory内存验收链；断言失败返回非零。
demo默认在测试临时目录建库，实际接纳样例、发布、登记更正、重开并导出旧成果；
它用Fraction从JSON证据独立复算，然后输出POINTWISE_SCENARIO单行JSON。
随后demo再演示连通链：偏移2的旧修订拒绝平移(1,0,0)，加入偏移3键的新修订
接纳同一平移，关闭重开后读取不漂移，并在离库后仅用json、hashlib核验
两份导出证据，输出第二条POINTWISE_SCENARIO。
不需要安装、联网或准备数据。运行清单见.pointwise-runtime.json。

对已有数据库演示指定修订：

```sh
python3.11 -B demo.py --database ./lab.sqlite --product specimen --revision 1
```

此方式不装载样例或更正来源，只执行指定修订并新增核算成果。
默认临时demo不修改项目样例。数据库初始化不灌样例、不重置用户数据；
兼容V1库重复打开保持所有修订、当前指针和成果，不兼容模式抛ValueError。

## 公开API和输入

```python
import json
from crystal import Archive, Tasks, Results

with Archive("lab.sqlite") as archive:
    with open("examples/day_cell.json", encoding="utf-8") as stream:
        cell = json.load(stream)
    revision = archive.put("specimen", cell)
    result_id = Tasks(archive).run_geometry("specimen", revision)
    result = Results(archive).read(result_id)
    evidence_json = Results(archive).export(result_id)
    assert result["bonds"][0]["length_squared"].numerator == 61
```

Archive(path=":memory:")支持上下文管理和幂等close()。
put(product_id,cell)提交完整替换版本，返回从1起递增的产品内修订号；
相同内容再次提交也登记新修订，摘要不用于自动去重。
read_revision(product_id,revision)返回规范nm完整cell字典；
current_revision(product_id)只查询当前指针。
Tasks(archive).run_geometry(product_id,revision)返回已提交成果整数ID；
每次正常运行可新增成果，必须显式指定修订。
Results(archive).read(result_id)读取持久成果；
export(result_id)返回完整JSON文本，不是文件路径。

cell恰含unit、lattice、sites、bonds四个必填字段：

- unit为nm、angstrom或pm，长度比例为1、1/10、1/1000，规范单位为nm。
- lattice是按行书写的3×3数组，数学上列为基向量A，行列式非零；
  负行列式合法。只缩放晶格，不缩放分数坐标或周期偏移。
- sites是{id,fractional}列表。id为非空字符串，fractional为[0,1)内三元组。
- bonds是{id,u,v,offset}列表。端点引用本修订原子，offset严格为三个Python整数。
  键(id,u,v,t)连接u和整数平移t处的v。

有理数接受非布尔整数、Fraction或两个非布尔整数构成的列表/元组[n,d]，
d不得为0。浮点、bool、None、缺字段及多余字段拒绝。
偏移不接受Fraction或有理数对，即使其值为整数。
产品、原子、键ID按原字符串比较，不裁剪或大小写折叠；原子和键分别要求唯一。
原子列表与键列表可为空，但任何键端点必须存在。零位移键合法。
同坐标不同原子、重复端点不同键ID、自像键、反向键均独立保留。
数值输出按身份排序。规范cell及所有JSON中的有理数为约分正分母整数对；
周期偏移、版本和身份编号保持整数。

## 结果、口径和证据

read返回独立嵌套容器，数值是Fraction：
source={product_id,revision,digest,calculation}，calculation为primitive-v1；
units={length:nm,length_squared:nm2,volume:nm3}；
volume为|det A|；
positions逐项为{id,cartesian}，位置为A f；
bonds逐项为{id,u,v,offset,displacement,length_squared}，
位移为A(f_v+t-f_u)，长度平方为位移内积。
不使用浮点容差、开方或坐标去重。

手算样例A=[[1,1/2,0],[0,1,0],[0,0,2]]，
u=(0,0,0)，v=(1/2,1/2,0)，t=(0,1,0)：
体积2，v位置(3/4,1/2,0)，位移(5/4,3/2,0)，长度平方61/16。
这些值实际经过任务发布、读取和导出，tests/test_geometry.py逐项核对。

export包格式为crystal-primitive-evidence、version=1。
source包含product_id、revision、digest、完整规范cell与原始送件raw_cell；
result包含与read相同的完整字段，但Fraction编码为[n,d]。
原始送件保留单位和送件顺序，Fraction转换为整数对，输入元组转换为JSON列表。
摘要为规范cell的SHA-256：字段固定、原子和键按ID排序，
json.dumps(ensure_ascii=True,sort_keys=True,separators=(",",":"),allow_nan=False)，
ASCII字节；摘要不包含任务ID、成果ID、时间或当前指针。
源码tests/test_export.py包含只依赖json、hashlib、Fraction的独立复算参照，
不会调用产品的矩阵、几何或摘要函数。导出前验证原始与规范来源一致，
并验证原子、键的逐项覆盖和数值。

## 周期像连通成果（connectivity-v1）

Tasks(archive).run_connectivity(product_id,revision,queries)返回持久成果整数ID；
产品、修订、成果ID的合法性和缺失语义沿用V1。仅本加工要求所引修订至多16原子、
32键，不限制接纳、几何加工或历史成果。queries为0～16项的list，每项恰含
u、v、translation；u、v须为该修订原子ID；translation为三个严格Python整数
构成的list或tuple，不接受bool、Fraction、浮点或有理数对，整数大小不另作限制。
空查询合法；空原胞仅可接收空查询；重复查询逐项保留，不合并、不排序。
非法查询、端点或规模超限抛ValueError，整次不新增任务或成果。

原键(id,u,v,t)在每个n∈Z³处连接(u,n)与(v,n+t)，允许正反向行走；反向仍用
同一原键ID，贡献相反偏移和负使用次数，不改V1有向位移与键身份语义。
逐项精确判断(u,0)到(v,translation)是否存在有限走向：在有限无向图连通分量内
取生成树，树路径给出基础偏移，剩余偏移须落在弦环偏移生成的整数格上；
可解性由Smith标准型精确整数消元判定，不以实数/有理数张成、浮点容差或
固定展开窗口代替整数可达性，也不求最短路。仅有一个原子u和键a=(u,u,(2,0,0))
时平移(1,0,0)不可达、(2,0,0)可达；加入独立键b=(u,u,(3,0,0))后(1,0,0)可达，
证据可为{"a":-1,"b":1}。

read对新类型返回恰含source、answers的独立字典；source恰含product_id、
revision、digest、calculation，calculation为connectivity-v1，digest仍是
完整规范cell的既有摘要。answers按输入顺序排列，每项恰含query、reachable、
coefficients；query为独立规范查询字典，translation统一为list。
不可达时coefficients为None；可达时为原键ID到非零严格整数的dict，
不记零项，没有净使用的键时用空字典。证据任选，不要求最短、唯一或最小系数。
系数是原键正向净使用次数，须满足：u、v在原胞有限无向图同一连通分量
（孤立原子自身也是分量）；每个出现的键ID都属于该分量；每个原子系数的
流入减流出等于它是v的指示量减去它是u的指示量，自像键在端点收支上相消；
所有键offset乘系数后的三维整数和恰为translation。

### 表示、原键证据及固定成果交接不变量

- 表示：连通任务与成果沿用tasks/results/result_bodies三表，与V1共用
  一个成果ID序列；tasks.calculation的CHECK放宽为
  ('primitive-v1','connectivity-v1')。打开V1库时在初始化事务内重建
  tasks表并保留全部行，user_version保持1；旧修订、当前指针、几何成果
  与成果ID不变，重建后foreign_key_check为空，重复打开不重复迁移。
- 成果体payload为{"version":1,"result":{source,answers}}，完整保存
  源修订绑定、规范查询、布尔结论与整数系数证据。Results.read只解码
  持久成果，不重新求解；Results.export只核验原始与规范来源一致性、
  摘要、头体绑定及正例证据，不调用连通求解器重新生成结论。
- 任务、成果头和成果体在同一事务原子发布；计算或发布异常原样传播，
  不留本次任务或部分成果，旧成果不变。更正来源、修改原查询或返回
  容器、关闭重开均不得改变旧成果；离库核账只需导出包，不打开数据库。
- export包恰含format、version、source、result：format为
  crystal-connectivity-evidence、version为整数1；source沿用V1导出的
  product_id、revision、digest、完整规范cell和raw_cell；result与
  新类型read返回的完整字典相同，布尔、整数和None用原生JSON表示。
  V1导出格式及既有值相等语义不变。

## 错误与失败边界

非法产品、结构、数值、单位和参数抛ValueError。
合法但不存在的产品/正整数修订/成果抛KeyError。
计算或数据库写入异常原样传播；本次修订与指针、任务与成果各自原子发布，
失败不留下空产品或半成品。返回值、原始输入的后续修改不写回归档。
旧修订及成果不因更正、读取或重开改变。证据缺失或来源不一致拒绝导出。
API不提供删除旧来源、隐式latest、并发调度或崩溃恢复协议。
archive.connection是底层SQLite连接；外部SQL破坏数据不属于正常写入API，
但读取会拒绝缺失成果体和头体错配。测试用真实SQLite触发器在写入中途拒绝，
验证事务回滚，并不替换计算函数。

## 代码与回归

接纳：intake/input_model/exact_numbers/units/input_digest。
几何：lattice/rational_matrix/sites/cell_geometry。
周期图及核账：bonds/periodic_graph/bond_geometry/bond_audit。
修订：archive/revision_store/revision_codec/schema/database。
成果：tasks/executor/results/export/result_store/result_codec。
连通：connectivity_queries/connectivity_graph/connectivity_solver/
smith_normal_form/connectivity_result/connectivity_audit。
所有30个模块均接入公开API，输入和结果值对象使用不可变元组与冻结记录。

N1、N3在test_intake；N2在test_geometry；N4、N5在test_bonds；
N6在test_tasks；N7在test_revisions；N8在test_export。
test_codec补充相同业务边界的精度、完整字段和模式版本验证。
test_connectivity覆盖连通链：手算样例与独立证书核账、耦合整数格与大数
Bezout、身份/方向、输入与规模整批拒绝、发布失败回滚、成果固定与重开、
离库证据核账、V1库迁移接纳新任务、连通codec边界及Smith标准型性质。
完整套件还覆盖固定seed=1729的有限三角晶格、左右手晶格、
真实发布/更正中途失败、大整数有理数、重复初始化及修改返回容器后重开。
测试创建的数据库仅位于测试临时目录，结束后清理。
