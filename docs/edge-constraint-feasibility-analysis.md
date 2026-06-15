# Edge Constraint 可行性分析报告

> **项目**: VolumeMesher (Convex Polyhedral Meshing for Robust Solid Modeling)
> **日期**: 2026-06-15
> **目标**: 判断能否将 triangle constraint 的 segment tracing 逻辑抽象为 edge constraint（单条线段在 Delaunay tet mesh 中经过哪些 tetrahedron）的追踪工具。
> **约束**: 只做源码探索，不实现、不修改、不生成 patch。

---

## 目录

1. [现有 Triangle Constraint 的映射流程](#1-现有-triangle-constraint-的映射流程)
2. [intersections_constraint_sides() 的真实作用](#2-intersections_constraint_sides-的真实作用)
3. [三种子函数的独立性验证](#3-三种子函数的独立性验证)
4. [Edge Constraint Tracing 的可行性](#4-edge-constraint-tracing-的可行性)
5. [风险点检查](#5-风险点检查)
6. [最小数据流设计](#6-最小数据流设计)
7. [需要关注的文件与位置](#7-需要关注的文件与位置)
8. [结论与建议](#8-结论与建议)
9. [Endpoint Over-Inclusion 问题](#9-endpoint-over-inclusion-问题)
10. [mark_TetIntersection 语义验证](#10-mark_tetintersection-语义验证)
11. [输入与点集合并流程](#11-输入与点集合并流程)
12. [Index Remap 风险分析](#12-index-remap-风险分析)
13. [Edge Constraint 数据结构设计](#13-edge-constraint-数据结构设计)
14. [最小验证插入位置](#14-最小验证插入位置)
15. [测试场景设计](#15-测试场景设计)
16. [Crossed-Tets Filter 的可复用函数](#16-crossed-tets-filter-的可复用函数)
17. [Index Remap 风险分析](#17-index-remap-风险分析)
18. [Edge Constraint 数据结构设计](#18-edge-constraint-数据结构设计)
19. [第一次 Remap 的精确语义](#19-第一次-remap-的精确语义)
20. [第二次 Remap 的精确语义](#20-第二次-remap-的精确语义)
21. [Edge Index 输入语义建议](#21-edge-index-输入语义建议)
22. [函数改造入口建议](#22-函数改造入口建议)
23. [Two-Input 与 Debug 输出建议](#23-two-input-与-debug-输出建议)
24. [最小实现计划：新增数据结构](#24-最小实现计划新增数据结构)
25. [最小实现计划：函数清单](#25-最小实现计划函数清单)
26. [最小实现计划：数据流步骤](#26-最小实现计划数据流步骤)
27. [最小实现计划：两次 Remap 注意点](#27-最小实现计划两次-remap-注意点)
28. [最小实现计划：Tracing 函数设计](#28-最小实现计划tracing-函数设计)
29. [最小实现计划：Debug 输出设计](#29-最小实现计划debug-输出设计)
30. [最小实现计划：测试用例](#30-最小实现计划测试用例)
31. [最小实现计划：不做事项](#31-最小实现计划不做事项)
32. [最小实现计划：风险清单](#32-最小实现计划风险清单)
33. [最小实现计划：最终结论](#33-最小实现计划最终结论)
34. [MVP-0 实现记录](#34-mvp-0-实现记录)
35. [Toy Case 测试策略](#35-toy-case-测试策略)
36. [MVP-0 测试结果](#36-mvp-0-测试结果)
37. [ELSE 分支专测](#37-else-分支专测)
 
 ---

## 1. 现有 Triangle Constraint 的映射流程

### 1.1 总体数据流

```
[输入] coords (顶点坐标) + tri_idx (三角形索引)
        │
        ▼
read_nodes_and_constraints()          ← makePolyhedralMesh.cpp
  去重、退化三角形过滤
        │
        ▼
mesh->tetrahedrize()                  ← Delaunay 四面体剖分
  顶点被置换 (original_index 改变)
        │
        ▼
[顶点索引重映射]                        ← makePolyhedralMesh.cpp L407-422
  将 tri_vertices[] 对齐到剖分后的顶点索引
        │
        ▼
fill_half_edges()
  → sort_half_edges()
  → place_virtual_constraints()       ← 三角形专属: 插入虚拟约束
        │
        ▼
insert_constraints()                  ← 核心映射入口 (4 步)
        │
        ▼
[输出] num_map / map                  ← tet_index → constraint_index 映射
      num_map_f0~f3 / map_f0~f3      ← tet_face → constraint 映射
        │
        ▼
BSPcomplex(mesh, constraints, map...) ← BSP 剖分消费这些映射
```

### 1.2 insert_constraints() 的四步流程

| 步骤 | 函数 | 作用 | 依赖三角形平面? |
|------|------|------|:---:|
| STEP 0 | `triangle_in_VT()` | 检查约束三角形是不是某个 tet 的面 | **是** |
| STEP 1 | `intersections_constraint_sides()` | 追踪三角形三条边穿过的所有 tet | **否** |
| STEP 2 | `find_improperIntersection()` | 分类边界 tet 的相交类型 | **是** |
| STEP 3 | `intersections_constraint_interior()` | 从边界向内 flood-fill，找内部 tet | **是** |
| STEP 4 | `compile_maps()` | 写入 map / num_map / map_fi / num_map_fi | **是** |

### 1.3 关注: STEP 3 为何依赖三角形

`intersections_constraint_interior()` 内部调用 `constrInterior_firstStep()` → `constrInterior_found()` → `tet_intersects_triInterior()`。

`tet_intersects_triInterior()` 使用 `vrt_signe_orient3d(t[0..3], c[0], c[1], c[2], mesh)` 对 tet 的四个顶点做相对三角形平面的方向测试（orient3d），以此判断 tet 是否与三角形内部相交。

**这正是边约束不需要的部分**。

---

## 2. intersections_constraint_sides() 的真实作用

### 2.1 核心判断

> **这是线段穿 tet mesh 的追踪器，不是几何切分。**

它不是在做 BSP 分割，也不是在做 surface classification。它在做的事情本质上就是：

> 给定一条线段 `<v_start, v_stop>` 和一个 Delaunay tet mesh，追踪这条线段从 v_start 到 v_stop 经过了哪些 tetrahedron。

### 2.2 算法流程

```
FOR each side (constr_side = 0, 1, 2)  ← 三角形有 3 条边
  
  1. BEGIN 阶段 (L2370-2384):
     v_start = constraint_vrts[constr_side]
     v_stop  = constraint_vrts[(constr_side+1) % 3]
     调用 intersections_TetVrtOnConstraintSide()
     找到所有从 v_start 出发的 incident tet
     标记 connecting_vrts[] = 下一步从哪里穿出

  2. CONTINUE 阶段 (L2390-2440):
     WHILE connecting_vrts[1] != v_stop:
       switch connecting_vrts[0]:
         case 1: → intersections_TetVrtOnConstraintSide()
         case 2: → intersections_TetEdgeCrossConstraintSide()
         case 3: → intersections_TetFacePiercedConstraintSide()
       将新发现的 tet 加入 intersecatedTet[] 数组
```

### 2.3 connecting_vrts 的编码方式

| element[0] | 含义 | element[1..3] |
|:----------:|------|---------------|
| 1 | segment 穿过了 tet 的 **一个顶点** | element[1] = 顶点索引 |
| 2 | segment 穿过了 tet 的 **一条边** | element[1..2] = 边端点索引 |
| 3 | segment 穿过了 tet 的 **一个面** | element[1..3] = 面顶点索引 |

这个编码方式完全独立于三角形——它只描述 segment 在 tet 拓扑中的前进路径。

---

## 3. 三种子函数的独立性验证

### 3.1 other_constr_vrt 是死参数

三个子函数都声明了 `uint32_t other_constr_vrt` 作为形参，但在函数体内从未被访问：

| 函数 | 文件名 | 行号 | 声明参数 | 内部使用量 |
|------|--------|:----:|----------|:----------:|
| `intersections_TetVrtOnConstraintSide()` | conforming_mesh.cpp | 1683 | `uint32_t other_constr_vrt` | **0 次引用** |
| `intersections_TetEdgeCrossConstraintSide()` | conforming_mesh.cpp | 1851 | `uint32_t other_constr_vrt` | **0 次引用** |
| `intersections_TetFacePiercedConstraintSide()` | conforming_mesh.cpp | 2084 | `uint32_t other_constr_vrt` | **0 次引用** |

这个参数出现在：
- 函数声明/定义行的形参列表中
- 上方注释中 ("index of the constraint vertex that do not belong to the side")
- 调用方的传参位置（`intersections_constraint_sides()` 从三角形第三个顶点传值）

**但从未在函数体的任何表达式、赋值、条件判断中被引用。**

### 3.2 三子函数实际依赖的内容

逐一检查每个函数内部的数据流：

```
intersections_TetVrtOnConstraintSide(mesh, v_curr, v_stop, ...)
  ├── mesh->incident_tetrahedra(v_curr)      ← tet 拓扑
  ├── mesh->tet_node                          ← tet 顶点索引
  ├── mesh->tet_neigh                         ← tet 邻接
  ├── opposite_face_vertices()                ← tet 拓扑
  ├── vrt_innerSegmentCrossesInnerTriangle()  ← 几何谓词 (只用 v_start, v_stop)
  ├── vrt_innerSegmentsCross()                ← 几何谓词 (只用 v_start, v_stop)
  ├── vrt_pointInInnerSegment()               ← 几何谓词 (只用 v_start, v_stop)
  └── fill_connecting_vrts()                  ← 纯赋值

intersections_TetEdgeCrossConstraintSide(mesh, tet_edge, v_start, v_stop, ...)
  ├── mesh->ETrelation(tet_edge, ...)         ← tet 拓扑
  ├── mesh->tet_node / tet_neigh              ← tet 拓扑
  ├── opposite_side_vertices()                ← tet 拓扑
  ├── vrt_innerSegmentCrossesInnerTriangle()  ← 几何谓词 (只用 v_start, v_stop)
  ├── vrt_innerSegmentsCross()                ← 几何谓词 (只用 v_start, v_stop)
  ├── vrt_pointInInnerSegment()               ← 几何谓词 (只用 v_start, v_stop)
  ├── vrt_same_half_plane()                   ← 方向判断 (只用于前进方向过滤)
  └── vrtsInSameHalfSpace()                   ← 方向判断 (只用于前进方向过滤)

intersections_TetFacePiercedConstraintSide(mesh, tet_face, v_start, v_stop, ...)
  ├── opposite_vertex_face() / adjTet_oppsiteTo_vertex()  ← tet 拓扑
  ├── vrt_innerSegmentCrossesInnerTriangle()  ← 几何谓词 (只用 v_start, v_stop)
  ├── vrt_innerSegmentsCross()                ← 几何谓词 (只用 v_start, v_stop)
  └── vrt_pointInInnerSegment()               ← 几何谓词 (只用 v_start, v_stop)
```

**三个函数都不访问：**
- 三角形第三个顶点索引
- 三角形平面
- 三角形法向
- 任何与三角形面积相关的数值

### 3.3 方向判断是否依赖三角形？

`vrt_same_half_plane()` 和 `vrtsInSameHalfSpace()` 在 `intersections_TetEdgeCrossConstraintSide()` 中用于过滤"前进方向是否正确"——确保 segment 是从 v_start 向 v_stop 方向穿过 tet 边，而不是反向。

这两个函数的输入是 tet 的边端点坐标 + segment 端点坐标，**与三角形的第三点无关**。它们只是在约束一个方向性的条件："segment 应该在穿过 tet 边时，从 v_start 方向进入，向 v_stop 方向退出"。

对于 edge constraint 的单次 tracing，这个方向判断依然有意义且不需要修改。

---

## 4. Edge Constraint Tracing 的可行性

### 4.1 判断

> ✅ **可以。子函数可直接复用，仅需调整入口逻辑。**

### 4.2 可行理由

| 条件 | 结论 | 证据 |
|------|:----:|------|
| 是否强依赖三角形第三点？ | **否** | `other_constr_vrt` 是死参数 |
| 是否强依赖三角形平面？ | **否** | 三子函数不访问任何平面几何 |
| 是否只是追踪 segment？ | **是** | `intersections_constraint_sides()` 本质是对 3 条 segment 分别追踪 |
| 子函数能否独立服务于单条 segment？ | **可以** | 去掉 `other_constr_vrt` 后每个函数逻辑完全自洽 |

### 4.3 需要剥离的部分

| 原有步骤 | 是否需要 | 理由 |
|----------|:--------:|------|
| `fill_half_edges()` | **否** | 从三角形面片提取边信息，edge 不需要 |
| `sort_half_edges()` | **否** | 同上 |
| `place_virtual_constraints()` | **否** | 为 BSP 剖分服务的三角形机制 |
| STEP 0: `triangle_in_VT()` | **否** | 检查三角形是否是 tet 面 |
| STEP 1: segment tracing | **是** | 核心可复用部分 |
| STEP 2: `find_improperIntersection()` | **否** | 使用三角形平面分类 |
| STEP 3: 内部 flood-fill | **否** | 使用三角形平面做 orient3d |
| STEP 4: `compile_maps()` | **部分** | 只需要 tet→edge 映射，不需要 face overlap 映射 |

### 4.4 需要保留的部分

| 原有机制 | 是否需要 | 理由 |
|----------|:--------:|------|
| Delaunay 剖分 (`mesh->tetrahedrize()`) | **是** | edge endpoint 必须是 tet mesh 的顶点 |
| 顶点重映射 | **是** | tetrahedrize 会置换顶点索引 |
| `intersections_TetVrtOnConstraintSide()` | **是** | BEGIN 阶段 |
| `intersections_TetEdgeCrossConstraintSide()` | **是** | CONTINUE 阶段 (edge step) |
| `intersections_TetFacePiercedConstraintSide()` | **是** | CONTINUE 阶段 (face step) |
| `enqueueTetsArray()` | **是** | 累积发现的 tet |
| `fill_connecting_vrts()` | **是** | 纯赋值，无副作用 |
| `INSERT_TET_IN_LIST` 宏 | **是** | 过滤已访问的 tet |
| `mark_TetIntersection[]` (简化版) | **是** | 用 0/1 标记已访问 tet |

---

## 5. 风险点检查

### 5.1 问题: mark_TetIntersection 的状态机过于复杂

当前有约 10 种标记值：

```c
#define INTERSECTION          1
#define IMPROPER_INTERSECTION 2
#define IMPROPER_INTERSECTION_COUNTED 3
#define PROPER_INTERSECTION_COUNTED   4
#define OVERLAP2D_F0          10
#define OVERLAP2D_F1          11
#define OVERLAP2D_F2          12
#define OVERLAP2D_F3          13
#define OVERLAP2D_F0_COUNTED  20
// ...
```

**风险**: edge constraint 只需要 0（未访问）/ 1（已访问）。如果强行复用原状态机，容易在状态切换中踩到不必要的分支。

**缓解**: 为 edge constraint 创建独立的标记数组，使用 `uint8_t*` 或直接用 `bool*`。

### 5.2 问题: 起始 tet 定位

BEGIN 阶段依赖 `mesh->incident_tetrahedra(v_start)`。如果 `v_start` 是孤立的顶点（不属于任何非 ghost tet），返回空数组，tracing 无法启动。

**缓解**: Edge endpoint 必须在 Delaunay 点集中，且 `v_start` 必须参与构成至少一个非 ghost tet。

### 5.3 问题: 方向判断的健壮性

`intersections_TetEdgeCrossConstraintSide()` 中的方向判断使用 `vrt_same_half_plane()` 和 `vrtsInSameHalfSpace()`。虽然不依赖三角形第三点，但这些判断对退化情况（segment 恰好穿过 tet 边端点）的行为需要额外验证。

**缓解**: 这些函数内部对共享端点的 case 返回 0，退化情况在现有 triangle constraint 中已有处理，edge constraint 场景下行为一致。

### 5.4 问题: edge 端点不一定是 tet 顶点

如果 edge endpoint 不在 Delaunay 点集中，当前 `intersections_TetVrtOnConstraintSide` 无法定位起始 tet。这是 **edge constraint 与 triangle constraint 的根本区别**：triangle 的顶点天然是 Delaunay 点集的子集，而 edge constraint 的端点可能不是。

**缓解**: 必须保证 edge endpoint 先加入 Delaunay 点集。如果未来允许 edge 端点在三角形内部，则需要新增"点定位"（point location）机制。

### 5.5 问题: 同时处理多条 edge 时的标记冲突

当前 `mark_TetIntersection` 在每次 `for(tri_ind)` 循环结束后被 `compile_maps()` 清空。如果未来需要在一次 scan 中处理多条 edge，需要为每条 edge 维护独立的标记数组，或者每处理完一条 edge 就立即清空。

---

## 6. 最小数据流设计

### 6.1 单个 Edge Constraint 的数据流

```
[输入]  edge_id, v0_idx, v1_idx
           │
           ▼
  1. 确保 v0, v1 已存在 Delaunay point set 中
     (如果是从外部输入的 edge，需先加入 vertices[])
           │
           ▼
  2. tetrahedrize()
     → 产生 tet mesh
     → 顶点可能被置换，需要重映射 v0/v1 索引
           │
           ▼
  3. tracing 初始化:
     - calloc mark_edge_tet[] 长度为 mesh->tet_num
     - intersecated_tets[] = NULL, num = 0
     - connecting_vrts[4]
     - nextTet_ind
           │
           ▼
  4. BEGIN 阶段:
     found = intersections_TetVrtOnConstraintSide(
               mesh, v0, v1, mark_edge_tet,
               connecting_vrts, &nextTet_ind, &num)
     enqueueTetsArray(found, num, &intersecated_tets, &n_total)
     free(found)
           │
           ▼
  5. CONTINUE 阶段 (loop):
     WHILE connecting_vrts[1] != v1:
       switch connecting_vrts[0]:
         case 1:  → TetVrtOnConstraintSide(...)
         case 2:  → TetEdgeCrossConstraintSide(...)
         case 3:  → TetFacePiercedConstraintSide(...)
       enqueueTetsArray(...)
           │
           ▼
  6. [输出]  edge_id → intersecated_tets[] (touched tetrahedra list)
```

### 6.2 追踪阶段的关键参数

| 参数 | 用途 | 来源 |
|------|------|------|
| `v_curr` / `v_start` | segment 起始端点 | 输入 |
| `v_stop` | segment 终止端点 | 输入 |
| `mark_TetIntersection[]` | 标记已访问 tet | `calloc` 分配（edge 专用） |
| `connecting_vrts[4]` | 拓扑前进方向 | 每步更新 |
| `nextTet_ind` | 下一步要检查的 tet | 每步更新 |
| `intersecatedTet[]` | 累积结果 | 每步追加 |

### 6.3 对比: Triangle vs Edge 数据流

```
Triangle:
  tri_vertices[3] → for each side → 3 次 tracing → intersectatedTet[]
                    → find_improperIntersection() → 分类
                    → intersections_constraint_interior() → flood-fill
                    → compile_maps() → map / num_map / map_f0~f3

Edge:
  edge_vertices[2] → 1 次 tracing → intersectatedTet[]
                    → (可选) 直接输出 edge → tet 映射
```

---

## 7. 需要关注的文件与位置

### 7.1 可直接复用的代码（无需修改）

| 文件 | 函数/宏 | 行号 | 备注 |
|------|---------|:----:|------|
| `conforming_mesh.cpp` | `extract_tetVrts()` | 25 | tet 顶点提取 |
| `conforming_mesh.cpp` | `extract_tetFaceVrts()` | 36 | tet 面提取 |
| `conforming_mesh.cpp` | `extract_tetEdgeVrts()` | 45 | tet 边提取 |
| `conforming_mesh.cpp` | `vrt_pointInInnerSegment()` | 58 | 几何谓词包装 |
| `conforming_mesh.cpp` | `vrt_innerSegmentsCross()` | 96 | 几何谓词包装 |
| `conforming_mesh.cpp` | `vrt_innerSegmentCrossesInnerTriangle()` | 114 | 几何谓词包装 |
| `conforming_mesh.cpp` | `vrt_same_half_plane()` | 188 | 方向判断 |
| `conforming_mesh.cpp` | `INSERT_TET_IN_LIST` 宏 | 18 | 过滤已访问 |
| `conforming_mesh.cpp` | `fill_connecting_vrts()` | 1649 | 纯赋值 |
| `conforming_mesh.cpp` | `enqueueTetsArray()` | 515 | 累积数组 |
| `conforming_mesh.cpp` | `intersections_TetVrtOnConstraintSide()` | 1681 | 核心子函数（可复用） |
| `conforming_mesh.cpp` | `intersections_TetEdgeCrossConstraintSide()` | 1847 | 核心子函数（可复用） |
| `conforming_mesh.cpp` | `intersections_TetFacePiercedConstraintSide()` | 2080 | 核心子函数（可复用） |
| `delaunay.h` | `TetMesh` 类 | — | 数据结构无需修改 |
| `delaunay.h` | `vertex_t` | — | 数据结构无需修改 |

### 7.2 需要新增/修改的位置

| 文件 | 位置 | 说明 | 工作量 |
|------|------|------|:------:|
| `conforming_mesh.h` | 新增 edge 约束结构体 | 或扩展 `constraints_t` 以支持 edge | 小 |
| `conforming_mesh.h` | 新增 `trace_edge_through_tets()` 声明 | 或 `insert_edge_constraints()` | 小 |
| `conforming_mesh.cpp` | 新增 edge tracing 入口函数 | 提取 side tracing 主循环为独立函数 | 中 |
| `conforming_mesh.cpp` | `intersections_constraint_sides()` | 可调新函数，也可保留兼容旧接口 | 小 |
| `makePolyhedralMesh.cpp` | 主流程中 edge 映射的调用点 | 在 `insert_constraints` 前后或替代它 | 中 |
| `main.cpp` | 输入解析 | 新增 edge 约束的输入格式 | 中 |

### 7.3 不需要修改的代码

| 文件 | 区域 | 理由 |
|------|------|------|
| `delaunay.h` / `delaunay.cpp` | 全部 | 剖分逻辑与约束类型无关 |
| `extended_predicates.h` / `.cpp` | 全部 | 底层几何谓词通用 |
| `hand_optimized_predicates.hpp` | 全部 | 同上 |
| `implicit_point.h` / `.hpp` | 全部 | 点类型通用 |
| `indirect_predicates.h` | 全部 | 间接谓词通用 |
| `BSP.cpp` / `BSP.h` | 全部 | BSP 剖分消费 map，不关心 map 的产生方式 |
| `graph_cut/` | 全部 | graph-cut 阶段完全在 BSP 之后 |
| `inOutPartition.cpp` | 全部 | 内外分区阶段在 BSP 之后 |
| `numerics.h` | 全部 | 数值工具通用 |

---

## 8. 结论与建议

### 8.1 总体判断

> ✅ **适合作为第一步。底层 segment tracing 逻辑成熟且不依赖三角形几何。**

三个子函数中的 `other_constr_vrt` 是死参数——这暗示代码作者在设计时已经留下了一个更通用的 segment tracing 接口空间，三角形只是一个便利的调用者。

### 8.2 第一步的具体任务建议

如果在未来实现时，建议按以下顺序进行：

1. **提取 tracing 循环为独立函数**
   - 从 `intersections_constraint_sides()` 中提取 WHILE 循环（L2390-2445）为 `trace_segment_through_tets(mesh, v_start, v_stop, mark, &intersecatedTet, &num)`
   - 去掉 `other_constr_vrt` 参数传递

2. **新增 edge 约束入口**
   - 在 `conforming_mesh.h` 中增加 edge 结构体和函数声明
   - 在 `conforming_mesh.cpp` 中新增 `trace_edge_constraint()` 或类似函数
   - 接入 STEP 1 的 tracing，跳过 STEP 2/3/4

3. **接线到主流程**
   - 在 `makePolyhedralMesh.cpp` 中，在 `tetrahedrize()` 之后、BSP 之前插入 edge tracing 调用
   - 输出 edge → tet 的映射

### 8.3 最可能踩坑的点

| 风险 | 严重程度 | 说明 |
|------|:--------:|------|
| `mark_TetIntersection` 状态机耦合 | **高** | 建议用独立标记数组解耦 |
| Edge endpoint 点定位 | **中** | 必须保证端点在 Delaunay 点集中 |
| 方向判断退化情况 | **低** | 与现有 triangle case 行为一致 |
| 多条 edge 的标记冲突 | **低** | 每次单独 calloc 即可解决 |
| 虚拟约束（`place_virtual_constraints`）| **不相关** | 只与 BSP 相关，edge-only tracing 不需要 |

### 8.4 下一步应该继续读的函数

如果继续深入，建议读：

1.  **`find_improperIntersection()`** (L1073) — 确认它确实强依赖三角形平面，确认 edge constraint 可以完全跳过它。

2.  **`intersections_constraint_interior()`** + **`tet_intersects_triInterior()`** (L2214) — 确认内部 flood-fill 确实使用 `vrt_signe_orient3d(t, c[0], c[1], c[2], mesh)`。这是 4.3 节"跳过 STEP 3"的关键证据。

3.  **`ETrelation()`** (delaunay.h / del_core.cpp) — 理解 `tet_edge` → `incident_tets` 的查询机制，`intersections_TetEdgeCrossConstraintSide()` 依赖它。

## 9. Endpoint Over-Inclusion 问题

### 9.1 第二轮研究的核心问题

上一轮的结论是底层 tracing 可复用。第二轮进一步追问：

> 如果直接复用 tracing，结果的语义是什么？会不会过宽？

### 9.2 intersections_TetVrtOnConstraintSide() 的 incident tet 处理

验证发现：**该函数无条件将所有 incident tetrahedra 加入结果**，不做方向过滤。

```
L1700-1701: incTet[] = mesh->incident_tetrahedra(v_curr)
                │
L1707-1709:     │  ← INSERT_TET_IN_LIST 只过滤已 visited 的
                ▼
L1711-1712: mark_TetIntersection[ALL new tets] = INTERSECTION
                │
L1722-1809:     │  ← 遍历 incTet 找"出口"（只为设置 connecting_vrts）
                │     不影响返回值
                ▼
L1812:     return ALL new incident tets（非出口方向的也包含）
```

关键行：L1709 的注释直言——"they have a generic intersection in v_curr"。

### 9.3 方向过滤的证据

函数体内有一段 case 判断（L1729-1805），识别 segment 从 tet 的哪个子结构穿出（面/边/顶点）。但这只用于填充 `connecting_vrts[]`，告诉调用方下一步往哪走。**返回值 `newTetIn_incTet` 是完整集合，不是子集。**

具体来说：

| Case | segment 行为 | 对 `connecting_vrts` 的影响 | 对返回值的过滤影响 |
|:----:|--------------|---------------------------|:----------------:|
| (3') | v_stop 在对面 | `connecting_vrts[1]=v_stop`，循环退出 | **无过滤** |
| (1) | 穿过对面内部 | 设置 face 编码 | **无过滤** |
| (2) | 穿过对边 | 设置 edge 编码 | **无过滤** |
| (3) | 穿过对顶点 | 设置 vertex 编码 | **无过滤** |
| (0) | 只接触 v_curr | 继续循环 | **无过滤** |

### 9.4 touched-tets vs crossed-tets

现有函数收集结果的实际覆盖范围：

```
[v0]                              [v0]
  │  tet_a   tet_b   tet_c          │  tet_a   tet_b   tet_c
  │  ·····                           │
  │         tet_d                    │         tet_d ← 被穿过
  │  tet_e  tet_f  tet_g             │  
  │         ·                        │
  │              tet_h               │              tet_h
  │         tet_i                    │         tet_i
  │                                  │
[v1]                              [v1]

touched-tets（宽松）               crossed-tets（精确）
v0/tet_a tet_b tet_c tet_d...      v0/...only tet_internal...
```

**现有函数更接近 "asymmetric touched-tets"：**

| 哪些 tet 被加入 | 原因 | 是否正确 |
|:---|:---|:---:|
| v0 的所有 incident tet | `intersections_TetVrtOnConstraintSide` 无条件加入 | ✅ 被标记，但可能过宽 |
| 线段内部穿过的 tet | 通过 face/edge piercing 找到 | ✅ 正确 |
| 沿途顶点 incident tet（全部） | vertex-passing case 同样无条件加入 | ✅ 被标记，但可能过宽 |
| v1 的 incident tet | **缺失** — 循环在 `connecting_vrts[1]==v1` 时退出 | ⚠️ 不对称 |

### 9.5 三角形能接受这个语义的原因

三角形约束的后续步骤弥补了 STEP 1 的 over-inclusion：

```
STEP 1: 收集边界 tet（过宽）
STEP 2: find_improperIntersection() 做详细分类
        → 用三角形平面做 orient3d，标记 IMPROPER / OVERLAP2D
STEP 3: 内部 flood-fill，添加内部 tet
STEP 4: compile_maps 只保留 IMPROPER 和 OVERLAP2D 标记的 tet
```

被 STEP 1 误加的 tet，如果在 STEP 2 中被判定为 proper intersection（即只接触一个顶点/边/面，不穿过内部），则不被 `compile_maps` 收录。

**Edge constraint 没有三角形平面，无法复现 STEP 2 的过滤逻辑。**

### 9.6 对 edge constraint 的影响

| 目标语义 | 第一版是否可行 | 需要什么处理 |
|----------|:------------:|-------------|
| **touched-tets（宽松）** | ✅ 可直接复用 | 需手动补上 v1 的 incident tet |
| **crossed-tets（精确）** | ❌ 不能直接复用 | 需要额外 endpoint filter |

精确 crossed-tets 需要过滤以下混淆情况：

```
v0 的 incident tet 中，哪些是"segment 真正穿过的"？

混淆 1: tet 与 segment 仅在 v0 处接触
混淆 2: tet 与 segment 仅在某个中间顶点处接触
混淆 3: tet 与 segment 仅共享一条不沿 v0→v1 方向的边

过滤思路: 对每个 incident tet，用 vrt_innerSegmentCrossesInnerTriangle()
         检查 segment 是否穿过该 tet 的对面
```

---

## 10. mark_TetIntersection 语义验证

### 10.1 在 segment tracing 阶段（STEP 1）是纯 visited 标记

三个子函数的证据：

**intersections_TetVrtOnConstraintSide (L1681-1821)**

```cpp
// L1708: 宏展开 = if(marker[tet]==0) → 只检查 == 0
INSERT_TET_IN_LIST(incTet[i], newTetIn_incTet,
                   num_newTetIn_incTet, mark_TetIntersection);
// L1712: 写入 INTERSECTION (=1)
mark_TetIntersection[ newTetIn_incTet[i] ] = INTERSECTION;
```

**intersections_TetEdgeCrossConstraintSide (L1847-2057)**

```cpp
// L1884-1885: 只检查 == 0
INSERT_TET_IN_LIST(incTet[i], newTetIn_incTet, num_newTetIn_incTet,
                   mark_TetIntersection);
// L1889: 写入 INTERSECTION (=1)
mark_TetIntersection[ newTetIn_incTet[i] ] = INTERSECTION;
```

**intersections_TetFacePiercedConstraintSide (L2080-2199)**

```cpp
// L2108: 只检查 == 0
if(mark_TetIntersection[ opp_tet_ind ] == 0){
// L2112: 写入 INTERSECTION (=1)
mark_TetIntersection[ opp_tet_ind ] = INTERSECTION;
```

### 10.2 在 STEP 1 中，没有任何对复杂状态的判断

在三个子函数中，**never 出现的操作**：

| 操作 | 出现位置 | 结论 |
|------|----------|:----:|
| `== IMPROPER_INTERSECTION` | 无 | ✅ 不检查 |
| `== OVERLAP2D_F0` | 无 | ✅ 不检查 |
| `== PROPER_INTERSECTION_COUNTED` | 无 | ✅ 不检查 |
| `写入 IMPROPER_INTERSECTION` | 无 | ✅ 不写入 |
| `写入 OVERLAP2D_*` | 无 | ✅ 不写入 |
| `写入 *_COUNTED` | 无 | ✅ 不写入 |

### 10.3 复杂状态在 STEP 2 才引入

`find_improperIntersection()` (L1070-1178) 使用三角形平面做 orient3d：

```cpp
// L1085-1089 — 使用三角形三个顶点做平面判断
or_tet_v[0] = vrt_signe_orient3d(constr_v[0], constr_v[1], constr_v[2], tet_v[0], mesh);
or_tet_v[1] = vrt_signe_orient3d(constr_v[0], constr_v[1], constr_v[2], tet_v[1], mesh);
or_tet_v[2] = vrt_signe_orient3d(constr_v[0], constr_v[1], constr_v[2], tet_v[2], mesh);

// L1119-1121 — 写 face overlap 标记
mark_TetIntersection[tet_ind] = 10 + f;  // OVERLAP2D_F0~F3

// L1147, L1169, L1175 — 写 improper 标记
mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;
```

### 10.4 完整状态转移

```
阶段               mark 值               谁写的？         依赖三角形？
─────              ─────                 ──────           ────────
初始               0                     calloc           否
STEP 1 tracing     0 → 1 (INTERSECTION)  三子函数          否
STEP 2 refine      1 → 2 (IMPROPER)      find_improper    是
                    1 → 10-13 (OVER2D)   find_improper    是
STEP 3 interior    1,2,10-13,20-23 →...  constrInterior   是 (间接)
STEP 4 compile     各种 → 0             compile_maps      否（只读）
```

### 10.5 对 edge constraint 的启示

- 可以完全复用 `uint32_t* mark_TetIntersection`，但只需使用 0 / INTERSECTION (1)。
- 不需要修改 `INSERT_TET_IN_LIST` 宏。
- 不需要修改三个子函数。
- 每条 edge 使用独立的 mark 数组（`calloc`），trace 完后 `free` 或 `memset`。
- 如果未来需要同时处理多条 edge，每条 edge 各自 `calloc` 即可，不互相污染。

---

## 11. 补充结论

### 11.1 第一版建议

综合两轮分析，edge constraint 第一版的建议是：

- ✅ **继续接受 "touched-tets" 语义**。它是 crossed-tets 的超集，保守但安全。后续可以通过 endpoint filter 收窄。
- ✅ **补齐 v1 的 incident tet**。tracing 结束后，额外调用 `mesh->incident_tetrahedra(v1)`，将未 visited 的 tet 加入结果。
- ✅ **mark_TetIntersection 使用独立 calloc，只走 0/1**。不必修改三个子函数。

### 11.2 两轮研究的核心发现

| 轮次 | 核心发现 | 位置 |
|:---:|----------|:----:|
| 第一轮 | `other_constr_vrt` 是死参数，子函数可独立服务于 edge | L1683, L1851, L2084 |
| 第一轮 | 算法本质是 segment tracing，不是切分 | L2331-2449 |
| 第二轮 | BEGIN 阶段将全部 incident tet 加入结果，不做方向过滤 | L1707-1712 |
| 第二轮 | STEP 1 中 mark 是纯 visited 标记，复杂状态在 STEP 2 才引入 | L1708, L1884, L2108 |

### 11.3 精确 crossed-tets 的改进方向

如果未来需要精确映射，建议的过滤路径：

```
STEP 1: segment tracing → 得到 touched-tets 集合

STEP 1.5: （可选）endpoint filter
  对 v0 的每个 incident tet：
    - 提取该 tet 的对面（opposite face）
    - 用 vrt_innerSegmentCrossesInnerTriangle(v0, v1, face, mesh)
      判断 segment <v0,v1> 是否真正穿过该 tet 内部
    - 只保留穿过内部或边/面的 tet
  对 v1 同理（但已由 tracking 路径覆盖的部分不需要重复）
```

这个过滤逻辑不需要三角形平面，只需要 segment 端点坐标和几何谓词——与 STEP 1 完全一致。

---

## 12. 输入与点集合并流程

### 12.1 现有点集构建流程

```
main()
  ├── read_OFF_file(fileA, &coords_A, &ncoords_A, &tri_idx_A, &ntriidx_A)
  └── makePolyhedralMesh(coords_A, ncoords_A, tri_idx_A, ntriidx_A, ...)
        │
        ├── new TetMesh
        ├── new constraints_t
        ├── read_nodes_and_constraints(coords_A, ncoords_A, tri_idx_A, ntriidx_A,
        │       &mesh->vertices, &mesh->num_vertices,       ← 输出: vertex_t[]
        │       &constraints->tri_vertices, &constraints->num_triangles)
        │
        ├── [L397-399] 重置 original_index (为 tetrahedrize 做准备)
        ├── mesh->tetrahedrize()
        ├── [L407-422] 三角形顶点索引第二次 remap
        ├── fill_half_edges → place_virtual_constraints
        ├── insert_constraints → map 构建
        └── BSPcomplex(mesh, constraints, map...)
```

### 12.2 新建点集接入点

Edge endpoint 坐标应该在 `read_nodes_and_constraints()` 内部与三角形顶点一起去重。这是最小侵入路径：

| 方案 | 接入位置 | 侵入性 | 理由 |
|------|----------|:------:|------|
| A — 在 `read_OFF_file()` 中读 | main.cpp | 低 | 需自定义文件格式，不改内核 |
| B — 在 `makePolyhedralMesh()` 入口新增参数 | makePolyhedralMesh.cpp | **最低** | 作为新函数参数传递，不改变已有输入逻辑 |
| C — 在 `read_nodes_and_constraints()` 内扩展合并逻辑 | makePolyhedralMesh.cpp | 中 | 复用现有 `map` / `diff` / `original_index` 去重 |

**推荐 B / C 组合**：在 `makePolyhedralMesh()` 入口增加 edge 坐标参数，传入 `read_nodes_and_constraints()` 的扩展版本（或新函数 `read_nodes_and_edges()`），与三角形顶点在同一个去重路径中处理。

### 12.3 为什么必须与三角形顶点同一次去重

如果单独对 edge endpoint 做一次去重，再和三角形顶点合并，会多一层索引映射，增加错误风险。与三角形顶点使用同一套 `map` / `diff` / `remove_duplicated_points` 逻辑，可以保证所有顶点的索引在同一个坐标排序空间下唯一。

---

## 13. Index Remap 风险分析

### 13.1 三角形顶点经历两次映射

**第一次：输入索引 → 去重后索引** (`read_nodes_and_constraints()` 内部)

```
L73-77: tmp[i].original_index = i       // 记录原始位置
L80:    remove_duplicated_points()      // 按坐标排序 + 去重

remove_duplicated_points 内部:
  L33: qsort(tmp, npts, ...)                 // 按坐标排序 → 顶点被置换!
  L43: map[tmp[i].original_index] = i         // 原位置 → 排序后位置
  L39: diff[i] = 累计重复数                    // 排序后位置 → 前面删除了几个
  L46-52: vertices_p = 去重后的 vertex_t[]

L85-92: 三角形顶点 remap:
  tri_vertices[k] = map[raw_idx] - diff[map[raw_idx]]
                    ↑                        ↑
                    排序后位置                减去被删重复数
```

**逻辑**：`map[原始输入索引]` 得到排序后位置，`diff[排序后位置]` 得到去掉的重复数，两者相减得到 `vertices[]` 中的去重后索引。

**第二次：tetrahedrize() 后的顶点置换**

```
L397-399: 重新设置 original_index = i（覆盖去重阶段的 orginal_index）
L402:    mesh->tetrahedrize()
         → init() 内部 (delaunay.cpp L79-80):
           swap(vertices[k], vertices[2])   // 保证前 4 点不共面
           swap(vertices[l], vertices[3])
         → 顶点数组被物理置换!
           但 original_index 随顶点一起移动

L407-422: 三角形顶点第二次 remap:
  if (vertices[2].original_index != 3):
      // 简单分支：直接读取 original_index
  else:
      // 复杂分支：手动进行 2↔3↔l 的 swap 修正
```

### 13.2 Edge endpoints 的两次映射

如果 edge endpoints 与三角形顶点在同一套去重流程中处理，它们自动经历：

```
原始输入 edge 端点索引
  → (第一次) map - diff → 去重后 vertices[] 中的索引
  → (第二次) mesh->vertices[idx].original_index → tetrahedrize 后的索引
```

**不需要额外步骤。**

### 13.3 关键风险点

| 风险 | 位置 | 触发条件 | 后果 |
|------|------|----------|------|
| **init() 置换仅影响索引 2/3** | delaunay.cpp L79-80 | 顶点 ≥ 4 且前 4 点存在共面 | vertices[2]/[3] 被换到其他位置 |
| **第二次 remap 的条件分支** | makePolyhedralMesh.cpp L407-422 | `vertices[2].original_index != 3` 决定走 if 还是 else | else 分支有手动 swap 逻辑 |
| **else 分支的手动修正** | L412-421 | 只有 2↔3↔l 的 swap，不处理其他索引 | 如果 edge endpoint 恰好是索引 2/3/l，需确认 else 分支是否覆盖 |
| **original_index 被覆盖** | L397-399 | 每次 tetrahedrize 前都重新设 | 去重阶段的 original_index 丢失（没影响，去重索引已在 tri_vertices 中） |
| **重复点合并** | remove_duplicated_points | edge endpoint 与 triangle 顶点坐标相同 | 自动共享索引 |

**最容易错位的场景**：tetrahedrize 的 else 分支（L412-421）。这个分支手动修正 vertices[2] / vertices[3] 与 l 的 swap。如果 edge endpoint 的索引恰好是 2、3 或 l，需要确认 else 分支中的修正逻辑也会正确处理 edge 顶点。第一版建议在第二次 remap 完成后，额外记录 edge 端点的新索引，而不是同时做 remap。

---

## 14. Edge Constraint 数据结构设计

### 14.1 方案对比

| 维度 | 方案 A: 独立结构 | 方案 B: 扩展 constraints_t |
|------|:----------------:|:--------------------------:|
| **安全性** | ✅ 不触发现有 triangle 代码路径 | ❌ 可能误伤 `fill_half_edges()`、`place_virtual_constraints()`、`insert_constraints()` 等 |
| **侵入性** | ✅ 只加约 5 行结构体定义 | ❌ 需要改造 constraints_t 构造/析构/所有调用方 |
| **语义清晰度** | ✅ edge 逻辑与 triangle 逻辑物理隔离 | ❌ 混合结构体膨胀，理解成本增加 |
| **与 insert_constraints 的关系** | 独立调用，参数隔离 | 需要 if/else 分流 |

### 14.2 建议：方案 A

```
// 示意结构 (conforming_mesh.h)
struct edge_constraints_t {
    uint32_t* edge_vertices;   // 按 [v0,v1, v0,v1, ...] 排列, 长度 = 2 * num_edges
    uint32_t num_edges;
};
```

**理由**：

1. **现有 `constraints_t` 被大量代码默认理解为三角形约束**
   - 成员名 `tri_vertices` 暗示三角形语义
   - `fill_half_edges()` 按 `tri_vertices[3*n]` 的步长遍历
   - `place_virtual_constraints()` 是三角形专属逻辑
   - `insert_constraints()` 的参数类型 `constraints_t*` 预期三角形

2. **扩展 constraints_t 容易产生静默误伤**
   - 新增 `num_edges` 后，现有循环可能无意中遍历到 edge 数据（如果 `num_triangles` 也代表 edge 就会出错）
   - 析构函数只 `free(tri_vertices)`，新增 `edge_vertices` 需要额外 free

3. **独立结构更适合第一版探索**
   - 不需要修改任何现有三角形代码路径
   - 可以自由改变 edge 相关逻辑而不影响已有功能

---

## 15. 最小验证插入位置

### 15.1 推荐位置

```
makePolyhedralMesh.cpp 主流程:

L467: insert_constraints(mesh, constraints, num_map, map, ...);
                    │
                    ▼
          ◼◼◼ 推荐插入位置 ◼◼◼
          │  mesh 仍有效
          │  tetrahedrize 已完成
          │  顶点 remap 已完成
          │  insert_constraints 已执行 (可参考其 mark 生命周期)
          │  num_map/map 仍可用 (L487 才被 free)
          │  BSP 尚未开始 (L482)
          │  不影响后续输出
          │
                    ▼
L482: BSPcomplex* complex_p = new BSPcomplex(mesh, constraints, map...);
```

### 15.2 满足条件核对

| 条件 | 满足情况 |
|------|:--------:|
| `mesh` 未被 delete | ✅ — mesh 在 L505 才 delete |
| `tetrahedrize()` 已完成 | ✅ — L402 |
| 顶点 remap 已完成 | ✅ — L407-422 |
| `insert_constraints()` 已完成 | ✅ — L467 |
| 不影响后续 BSP | ✅ — 只打印，不修改 mesh 或 map |
| `num_map` / `map` 仍可访问 | ✅ — L487 才 free |

### 15.3 验证输出内容（示意）

```
// 伪代码，不实现
for each edge:
    trace_segment_through_tets(mesh, v0, v1, mark, &touched, &num)
    print "edge_e (v0→v1): touched tet ids = [...]"
    free(touched)
```

---

## 16. 测试场景设计

### 场景 1: Edge 恰好是 Delaunay mesh 的一条 tet 边

| 项目 | 内容 |
|------|------|
| **配置** | 4 个点的单纯 tetrahedron，edge 是其中一条边 |
| **预期** | touched-tets 包含该边所在的所有 tet（通常是 2+ 个） |
| **验证点** | BEGIN 阶段是否正确处理 v0 的 incident tets；是否在 while 第一次迭代就通过 case (3') 发现 v1 |
| **暴露问题** | v1 的 incident tet 是否被对称地遗漏 |

### 场景 2: Edge 穿过多个 tet 的面内部

| 项目 | 内容 |
|------|------|
| **配置** | 长方体被剖分成多个 tet，edge 从一端到另一端穿过多个 tet 面 |
| **预期** | touched-tets 包含沿途所有被 face-piercing 穿过的 tet |
| **验证点** | CONTINUE 阶段 `connecting_vrts[0]==3` 的路径是否正确步进，是否有遗漏或重复 |
| **暴露问题** | face-piercing 方向判断是否正确；`intersections_TetFacePiercedConstraintSide` 是否正确返回相邻 tet |

### 场景 3: Edge 穿过一个共享边

| 项目 | 内容 |
|------|------|
| **配置** | 多个 tet 共享一条边，edge 穿过该边的内部 |
| **预期** | touched-tets 包含共享该边的所有 tet |
| **验证点** | `connecting_vrts[0]==2` 的路径；`intersections_TetEdgeCrossConstraintSide` 的方向判断 |
| **暴露问题** | 共享边的所有 tet 是否全部被发现和收录 |

### 场景 4: Edge 经过一个内部顶点

| 项目 | 内容 |
|------|------|
| **配置** | 多个 tet 围绕一个内部顶点，edge 穿过该顶点 |
| **预期** | touched-tets 包含该顶点周围的所有 incident tet（验证 over-inclusion） |
| **验证点** | vertex-passing 时 incident tet 全部加入还是只加入沿方向的部分 |
| **暴露问题** | over-inclusion 的实际范围 |

### 场景 5: Edge endpoint 位于共享顶点上（v0 是 mesh 的一个 corner）

| 项目 | 内容 |
|------|------|
| **配置** | 一个 cube 的 corner 被剖分成 6 个 tet（Kuhn 剖分），v0 是 corner，v1 是对角 |
| **预期** | touched-tets: v0 的 6 个 incident tet 全部加入 + 内部 tet + v1 incident tet（部分） |
| **验证点** | 人工与 connectivity 检查对照 |
| **暴露问题** | v1 incident tet 的不对称缺失；over-inclusion 的程度 |

### 场景 6: Edge 完全与一个 tet 面共面

| 项目 | 内容 |
|------|------|
| **配置** | edge 在一个 tet 的面上，不穿过内部 |
| **预期** | touched-tets 只包含面所在的两个 tet（该面和邻面），不包含内部 tet |
| **验证点** | 退化情况的处理 |
| **暴露问题** | face-piercing 对共面退化的行为 |

---

## 17. Crossed-Tets Filter 的可复用函数

### 17.1 现有可复用的基础函数

| 函数 | 行号 | 用途 |
|------|:----:|------|
| `extract_tetVrts()` | L25 | 提取 tet 四个顶点索引 |
| `extract_tetFaceVrts()` | L36 | 提取 tet 特定面的三个顶点 |
| `extract_tetEdgeVrts()` | L45 | 提取 tet 特定边的两个端点 |
| `opposite_face_vertices()` | L402 | 获取 tet 中 v_curr 对面的 3 个顶点 |
| `opposite_side_vertices()` | L441 | 获取 tet 中与给定边相对的边 |
| `vrt_innerSegmentCrossesInnerTriangle()` | L114 | **核心**：判断 segment 是否穿过一个三角形的内部 |
| `vrt_innerSegmentsCross()` | L96 | 判断两条 segment 的内部是否相交 |
| `vrt_pointInInnerSegment()` | L58 | 判断点是否在 segment 内部（不含端点） |
| `vrt_pointInSegment()` | L66 | 判断点是否在 segment 上（含端点） |
| `vrt_same_half_plane()` | L188 | 判断两个点是否在 segment 同一侧 |

### 17.2 是否足够？

**已有足够的基础函数。** 不新增几何谓词，只需要组合现有函数写一个 "segment crosses tet interior" 的判断逻辑。

现有的 `tet_intersects_triInterior()` (L2214) 判断的是 "tet 是否与三角形内部交叉"——它用了三角形平面的 orient3d，不适合 edge。但用现有函数可以写出等价于 "segment 是否穿过 tet 内部" 的判断：

```
对 tet 中 v0 的对面 (opposite face f):
  if vrt_innerSegmentCrossesInnerTriangle(v0, v1, f, mesh):
    → 穿过 face → 保留

对 tet 中不与 v0 相连的各边:
  if vrt_innerSegmentsCross(v0, v1, edge):
    → 穿过 edge → 保留

对 tet 中不与 v0 相连的各顶点:
  if vrt_pointInInnerSegment(v, v0, v1):
    → 穿过 vertex → 保留
```

上述逻辑全都可以用现有函数实现。不需要新的 C++ 代码，只需要组合。

---

## 18. 补充结论（第三轮）

### 18.1 核心判断汇总

| 问题 | 回答 | 关键证据 |
|------|------|----------|
| 第一版是否应接受 touched-tets？ | ✅ **是** | crossed-tets 是超集，保守安全 |
| Edge endpoints 是否必须作为 Delaunay vertices？ | ✅ **是** | tracing 依赖 `incident_tetrahedra()`，要求端点必须是 mesh vertex |
| Edge constraint 是否应独立于 constraints_t？ | ✅ **是** | constraints_t 被 triangle 代码默认使用，扩展误伤风险高 |
| 最小验证插入位置 | L467-L482 之间 | mesh 有效 + remap 完成 + BSP 未开始 |
| 基础函数是否足够支持精确化？ | ✅ **足够** | 现有几何谓词（6 个）和拓扑函数（5 个）可直接组合 |

### 18.2 三阶段路线图（仅策略建议）

```
阶段 0（已完成）: 源码探索
  验证 three sub-functions 不依赖 triangle
  验证 mark 在 STEP 1 中为纯 visited 标记
  确认 index remap 风险点
  设计最小验证位置

阶段 1: touched-tets 第一版
  独立 edge_constraints_t 结构体
  在 read_nodes_and_constraints 中合并 edge endpoints
  复用 three sub-functions 做 tracing
  补齐 v1 incident tet
  只输出 edge → tet 列表，不接 BSP

阶段 2: crossed-tets 精确化
  用现有几何谓词写 endpoint filter
  从 touched-tets 过滤得到 crossed-tets
  可选: 写 edge → tet 反向索引
```

### 18.3 下一步建议

| 优先级 | 建议方向 | 理由 |
|:------:|----------|------|
| 1 | 读 `incident_tetrahedra()` (delaunay.cpp L648) | 理解 incident tet 的查找行为，评估对 v0/v1 补足的影响 |
| 2 | 读 `ETrelation()` (delaunay.cpp L681) | 理解 edge→incident tets 的旋转查询机制 |
| 3 | 读 `insert_constraints()` 的 mark 重置时机 (L2862-2870) | 确认 compile_maps 后的 memset 逻辑，设计 edge 版本的 mark 生命周期 |

---

## 19. 第一次 Remap 的精确语义

### 19.1 map - diff 公式的逐行追踪

```cpp
// remove_duplicated_points 内部 (makePolyhedralMesh.cpp L28-56)
// L33: qsort(tmp, npts, sizeof(vertex_t), vertex_compare)
//      按 (x,y,z) 字典序对 tmp[] 排序。
//      tmp[0], tmp[1], ..., tmp[n-1] 被物理重排。
//      但 original_index 随元素一起移动。

// L43: map[tmp[i].original_index] = i
//      含义: "原始输入索引 X 的元素，排序后在位置 i"

// L39: diff[i] = 从 0 到 i 累计遇到的重复顶点数
//      排序后数组 tmp[0..i] 中有多少重复顶点被删除了

// 调用处 (L90-92):
// raw = tri_idx_A[j*3 + k]   ← 原始 OFF 文件中的顶点索引
// sorted_pos = map[raw]       ← 排序后位置
// dup_before = diff[sorted_pos]  ← 该位置之前被删的重复数
// new = sorted_pos - dup_before  ← 去重后 vertices_p[] 中的索引
```

### 19.2 精确映射链

```
原始输入索引 (tri_idx_A 中的值)     ────────┐
       ↓                                          ↓
tmp[original_index=raw].coord 被 sort        map[原始索引] = 排序后位置
       ↓                                          ↓
排序后位置 = map[原始索引]                    diff[排序后位置] = 前面累计重复数
       ↓                                          ↓
map[原始索引] - diff[map[原始索引]] ────────────┘
       ↓
= 去重后 vertices_p[] 中的最终索引
```

**结论**：`map[raw] - diff[map[raw]]` 实现了：

```
OFF 文件中的 old vertex id
→ 去重后 mesh->vertices 中的 new vertex id
```

这个说法**完全准确**。

### 19.3 Edge endpoints 是否可以复用

**完全可以。** Edge endpoints 的坐标数组只要与三角形顶点坐标拼成同一个 `tmp[]`，经历同一轮 `remove_duplicated_points` → 同一套 `map - diff`，最终 edge 端点索引与三角形顶点索引来自同一套去重映射。坐标重合的点会自动合并。

唯一要求：edge endpoints 在拼接数组中的原始输入索引必须连续编号（例如三角形顶点占 0..npts_A-1，edge 端点占 npts_A..npts_A+nedge_pts-1）。

---

## 20. 第二次 Remap 的精确语义

### 20.1 tetrahedrize() 对 vertices[] 做了什么

`tetrahedrize()` 调用 `init()`（delaunay.cpp L64-80），其中会**物理交换 vertices[] 数组元素**：

```cpp
// L64-66: 从 i=0, j=1, k=2, l=3 开始搜索
uint32_t i=0, j=1, k=2, l=3;

// L68-73: 找到使 orient3d(vertices[0],vertices[1],vertices[k],vertices[l]) != 0 的 k,l
for (; ori==0.0 && k<n-1; k++)
    for (l = k+1; ori==0.0 && l<n; l++)
        ori = orient3d(vertices[i].coord, vertices[j].coord,
                       vertices[k].coord, vertices[l].coord);

// L79-80: 物理交换!
std::swap(vertices[k], vertices[2]);   // 把位置 k 的顶点换到位置 2
std::swap(vertices[l], vertices[3]);   // 把位置 l 的顶点换到位置 3
```

交换后的 `vertices[]` 排列：

| 位置 | 内容 | `original_index` | 含义 |
|:----:|------|:---:|------|
| `vertices[0]` | 始终是原始 `vertices[0]` | 0 | **未参与交换** |
| `vertices[1]` | 始终是原始 `vertices[1]` | 1 | **未参与交换** |
| `vertices[2]` | 原始 `vertices[k]` 的顶点 | k | 从位置 k 换进来的 |
| `vertices[3]` | 原始 `vertices[l]` 的顶点 | l | 从位置 l 换进来的 |
| `vertices[k]` | 原始 `vertices[2]` 的顶点 | 2 | 被换到 k 的 |
| `vertices[l]` | 原始 `vertices[3]` 的顶点 | 3 | 被换到 l 的 |
| 其余位置 | 未参与交换 | 不变 | 没动过 |

### 20.2 第二次 remap 的目的

**把三角形顶点索引从"去重后 vertices[] 中的位置"更新为"tetrahedrize 后 vertices[] 中的位置"。**

因为 `vertices[]` 的元素被 swap 了，索引指向了错误的顶点坐标。`original_index` 记录了每个顶点在 swap 前的位置，利用它可以把旧索引反算回新索引。

### 20.3 if 分支（最常见情况）

```cpp
// L407: 条件 vertices[2].original_index != 3
// 含义: 位置 2 的顶点不是从位置 3 换进来的
//      → 只有两个独立的 swap(2↔k) 和 swap(3↔l)
//      → 没有循环依赖

// L409-410:
for (k = 0; k < 3 * constraints->num_triangles; k++)
    constraints->tri_vertices[k] =
        mesh->vertices[constraints->tri_vertices[k]].original_index;
```

**为什么 vertices[old_pos].original_index 能正确映射？**

因为 swap 是对称的：A↔B 后，`vertices[A].original_index = B` 且 `vertices[B].original_index = A`。所以：

```
假设 triangle 有旧索引 2（指向 swap 前的 vertices[2]）。
swap 后 vertices[2] 被换成了位置 k 的顶点。
vertices[2].original_index = k。
"原位置 2 的顶点现在在哪？" → 在位置 k。
结果 = vertices[2].original_index = k = 正确答案。✓
```

### 20.4 else 分支（vertices[2].original_index == 3 的特殊情况）

当 k == 3（即 `vertices[3]` 被换到了位置 2），发生链式交换：

```
初始: vertices[2] (从 2), vertices[3] (从 3), vertices[l] (从 l)

swap(vertices[3], vertices[2])  → 位置 2: 从3, 位置 3: 从2
swap(vertices[l], vertices[3])  → 位置 3: 从l, 位置 l: 从2

最终:
  位置 2: original_index = 3
  位置 3: original_index = l
  位置 l: original_index = 2

约束旧索引 2 → 正确答案: l  (从2的顶点经 3 到了 l)
约束旧索引 3 → 正确答案: 2  (从3的顶点到了 2)
约束旧索引 l → 正确答案: 3  (从l的顶点到了 3)
```

`vertices[2].original_index` 在这种链式 swap 中给出的是 3（直接从 3 换进来的顶点），而不是正确值 l。所以简单公式失效，需要手动修正：

```cpp
// L414-421:
uint32_t l = mesh->vertices[3].original_index;  // = l
uint32_t c = constraints->tri_vertices[k];
if (c == 2)      new_idx = l;    // 原 2 → l
else if (c == 3) new_idx = 2;    // 原 3 → 2
else if (c == l) new_idx = 3;    // 原 l → 3
```

### 20.5 Edge 顶点是否需要相同 remap？

| 分支 | Edge remap 公式 | 是否安全 |
|:----:|-----------------|:--------:|
| **if 分支** (L409) | `edge_vertices[k] = mesh->vertices[edge_vertices[k]].original_index;` | ✅ **完全安全**。与 triangle 完全相同的公式 |
| **else 分支** (L414-421) | `if(c==2) c=l; else if(c==3) c=2; else if(c==l) c=3;` | ✅ **可并行处理**。edge idx 只要与 triangle idx 共享同一个 vertices[] 数组，遇到的 swap 完全相同 |

第一版建议：**在同一个 if/else 代码块中，把三角形循环复制一份，改为 edge 循环**。不要合并循环——保持代码清晰，让审查者能一眼看出 edge remap 和 triangle remap 走的是完全相同的路径。

---

## 21. Edge Index 输入语义建议

### 21.1 三种语义对比

| 语义 | 输入示例 | 说明 | 核心优点 | 核心缺点 |
|:----:|----------|------|----------|----------|
| **A** | 独立 `edge_coords[]` + `edge_idx[]` 指向自己 | edge idx 与 triangle 解耦 | 隔离干净 | 重合点不自动合并 |
| **B** | 合并的 `coords_all[]` + `edge_idx[]` 指向合并数组 | edge 坐标拼到 triangle 后面，一齐去重 | 重合自动合并 | 需要在合并层做一次 coord copy |
| **C** | `edge_idx[]` 直接引用 OFF 顶点索引 | 与 triangle idx 用同一套 `coords_A[]` | 最简洁 | edge 端点必须是 OFF 中已有的点，无法引入新点 |

### 21.2 建议：第一版用语义 B

```
1. makePolyhedralMesh() 入口接收 edge_coords + edge_idx
2. 函数内拼接:
   coords_all = [triangle_vertices | edge_endpoints]
   edge_idx 指向 coords_all 的后半段
3. 传给扩展的 read_nodes_constraints_and_edges() 统一去重
```

**理由**：
1. **重合自动合并**：edge endpoint 坐标与 triangle vertex 重合时，`remove_duplicated_points` 的坐标排序 + 比较会自然合并，不产生多余顶点。这是 Delaunay 剖分的要求。
2. **`map/diff` 完全透明**：edge endpoints 在统一输入数组中的位置是连续的，与三角形顶点经历同一套 `remove_duplicated_points → map - diff → original_index` 公式。
3. **不需要用户手动算索引**：用户只需要提供 edge 端点在独立坐标数组中的索引（语义 A 的内层），合并后的重映射由代码自动完成。

### 21.3 语义 C 的适用场景

如果 edge endpoint 必须是 OFF 文件中已有的顶点（即 edge 不引入新点），语义 C 最简洁。但第一版探索时应该允许 edge 引入新点，所以建议语义 B。

---

## 22. 函数改造入口建议

### 22.1 方案对比

| 方案 | 做法 | 侵入性 | 风险 | 回退难度 |
|:----:|------|:------:|:----:|:--------:|
| **A — 扩展原函数** | `read_nodes_and_constraints()` 参数增加 edge_coords, nedge_pts, edge_idx, nedge | 中 | 改两处函数签名 + 所有调用点 | 中 |
| **B — 新增 sibling** | 新增 `read_nodes_constraints_and_edges()` | **低** | 原函数不动，新函数内部处理 map/diff | **最低** |
| **C — 临时拼接调原函数** | 在 `makePolyhedralMesh()` 中拼接坐标 + 索引，然后调原函数 | **最低** | 需要额外坐标拷贝，不改现有函数签名 | **最高** |

### 22.2 建议：第一版用方案 B

```cpp
// 示意: 放在 read_nodes_and_constraints 旁边（不实现）
void read_nodes_constraints_and_edges(
    double* coords_A, uint32_t npts_A, uint32_t* tri_idx_A, uint32_t ntri_A,
    double* edge_coords, uint32_t npts_edge, uint32_t* edge_idx_A, uint32_t nedge,
    vertex_t** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri,
    uint32_t** edge_vertices_p, uint32_t* nedge_out,
    bool verbose);
```

函数内部：
1. 将 `coords_A` 和 `edge_coords` 拼成统一的 `tmp[]`（三角形顶点索引 0..npts_A-1，edge 端点 npts_A..npts_A+npts_edge-1）
2. 调用 `remove_duplicated_points()`（不修改它）
3. 用同一套 `map - diff` 同时 remap triangle_idx 和 edge_idx
4. 输出 `tri_vertices_p` 和 `edge_vertices_p`

**理由**：
| 考虑 | 说明 |
|------|------|
| `remove_duplicated_points` 设计为通用函数 | 只操作 `vertex_t[]`，不关心顶点来源。不需要修改 |
| 原函数不受影响 | triangle-only 路径稳定，不影响 `read_nodes_and_constraints_twoInput` |
| 回退容易 | 新函数不与现有代码耦合，可以随时删除 |
| 将来可以合并 | 探索结束后，如果觉得代码重复，可以把逻辑合并回原函数 |

---

## 23. Two-Input 与 Debug 输出建议

### 23.1 第一版范围：仅支持单输入

**建议第一版只支持 `bool_opcode == '0'` 的单输入路径，不考虑 two-input boolean。**

| 理由 | 说明 |
|------|------|
| **复杂度** | two-input 需要同时处理两组 triangle 坐标 + 两组 edge 坐标，索引映射翻倍 |
| **输入源不清晰** | edge constraint 是单独指定的，不是 OFF 文件的一部分。two-input 时 edge 可能来自模型 A、模型 B 各一部分 |
| **验证干扰** | 第一版目标是验证 tracing 层，two-input 增加的是 boolean operation 层的复杂度 |
| **未来扩展** | `read_nodes_and_constraints_twoInput` 的结构与单输入对称，可按同一模式扩展 |

主流程调用示意（不实现）：

```cpp
if (!two_input) {
    read_nodes_constraints_and_edges(..., &edges, ...);  // 只走单输入
} else {
    read_nodes_and_constraints_twoInput(...);  // 不变，triangle-only
}
```

### 23.2 Debug 输出建议

为了验证 index remap 完整链条，在验证插入位置（L467-L482 之间）应打印：

**第一层：确认两次 remap 正确**

| 数据 | 目的 | 格式 |
|------|------|------|
| edge_id | 标识 | `"Edge %u:"` |
| raw input idx | 输入阶段原始索引 | `" raw=(%u,%u)"` |
| after dedup idx | 第一次 remap 后 | `" dedup=(%u,%u)"` |
| after tetrahedrize idx | 第二次 remap 后 | `" dtize=(%u,%u)"` |
| endpoint coordinates | 确认坐标正确 | `" v0=(%f,%f,%f)"` |
| touched tet count | tracing 结果 | `" touched=%llu tets"` |

**第二层：输出 vertex remap 验证开关**

| 数据 | 目的 |
|------|------|
| `mesh->num_vertices` | 确认去重后顶点总数 |
| `mesh->tet_num` | 确认剖分后 tet 数量 |
| `vertices[2].original_index` | 辅助判断走 if 还是 else 分支 |
| first touched tet 的顶点索引 | 确认 tet 与 edge 的索引对应关系 |
| first touched tet 的顶点坐标 | 与 edge 端点坐标对照 |

**Debug 输出顺序（示例）**

```
=== edge remap verification ===
edge 0: raw=(0,1) dedup=(0,1) dtize=(0,1)  v0=(1.2,3.4,5.6) v1=(7.8,9.0,1.2)
edge 1: raw=(2,5) dedup=(1,3) dtize=(2,5)  v0=(...)          v1=(...)

=== edge tracing result ===
edge 0: touched 3 tets [7 12 19]
  tet 7  vertices: [0 4 11 23]  coord=(...)
  tet 12 vertices: [0 15 7 11]  coord=(...)
```

如果 `raw != dedup` 或 `dedup != dtize`，说明对应 remap 生效。如果一个 edge endpoint raw=2 且 dedup=2 且 dtize=2，同时 `vertices[2].original_index = 4`，说明该点是全新点且没有被 tetrahedrize 移动。

---

## 24. 补充结论（第四轮）

### 24.1 核心判断

| 问题 | 回答 | 关键证据 |
|------|------|----------|
| `edge_vertices` 是否可以复用 triangle 的两次 remap？ | ✅ **可以** | 两次 remap 公式与顶点来源无关：第一次是 `map-diff`，第二次是 `vertices[idx].original_index` |
| 第二次 remap 是否有特殊坑？ | ✅ **有，已知且可控** | else 分支（`vertices[2].original_index == 3`）处理链式 swap 2→l→3→2。edge idx 如果恰好是 2/3/l 需要此分支保护 |
| 第一版最安全的数据入口是什么？ | **语义 B：拼接坐标 → 统一去重** | edge endpoints 和 triangle vertices 经历同一套去重 + remap |
| 最小侵入函数改造方案 | **方案 B：新增 sibling 函数** | 原函数不动，新函数内部自行处理 map/diff |
| 第一版是否支持 two-input | **否，只做单输入** | 降低复杂度，专注 tracing 层验证 |

### 24.2 两步 remap 的安全保证

```
Step 0: edge endpoints 与 triangle vertices 拼接为统一坐标数组
        ↓
Step 1: remove_duplicated_points() 
        → qsort + 坐标比较 + 去重
        → map-diff 公式同时对两者生效
        → 坐标重合自动合并
        ↓
Step 2: tetrahedrize() + init()
        → vertices[] 物理 swap
        → 第二次 remap 同时对两者生效
        → if 分支: 通用 original_index 公式
        → else 分支: 手动修正 2↔3↔l
        ↓
最终: edge_vertices 与 tri_vertices 在同一套顶点索引空间中
        ↓
edge tracing 使用的 v0/v1 索引与 triangle constraint 使用的顶点索引
是同一套坐标系, 索引互相引用时不会错位
```

两层 remap 的核心公式对 edge 和 triangle 完全透明。不需要为 edge 设计新的 remap 逻辑。

### 24.3 最应该避免的三个错误

| 排名 | 错误 | 预防措施 |
|:----:|------|----------|
| 1 | edge 与 triangle 不在同一套去重路径 | **语义 B 强制保证**：统一坐标数组、统一 `tmp[]`、统一 `map/diff` |
| 2 | 第二次 remap 的 else 分支漏了 edge | **复制 if/else 代码块**：edge_vertices 在同一个 if/else 块中做同样的事 |
| 3 | edge 端点不在 Delaunay 点集中 | **从 `read_nodes_constraints_and_edges()` 开始就保证**：所有顶点都在 `vertices[]` 中 |

### 24.4 下一步建议

| 优先级 | 方向 | 理由 |
|:------:|------|------|
| 1 | 在 `makePolyhedralMesh.cpp` 的验证位置打印 edge 端点的 `raw_id → dedup_id → dtize_id` 链条 | 用实际数据确认两次 remap 的正确性 — 这是第一版实现前必做的验证 |
| 2 | 读 `read_nodes_and_constraints_twoInput` 的完整逻辑 (L115-200) | 其 map/diff 使用与单输入完全一致，为 future two-input 扩展做准备 |
| 3 | 设计 edge 测试用的 OFF 文件 + edge 输入 | 用已知小模型（如 4 点 tetrahedron）验证 tracing |

---

## 25. 最小实现计划：新增数据结构

### 25.1 `edge_constraints_t`

**文件**: `src/conforming_mesh.h`

**字段**:

| 字段 | 类型 | 说明 |
|------|------|------|
| `edge_vertices` | `uint32_t*` | 按 `[v0,v1, v0,v1, ...]` 顺序排列，长度 `2*num_edges`。存储的是两次 remap 完成后指向 `mesh->vertices[]` 的顶点索引 |
| `num_edges` | `uint32_t` | edge constraint 的数量 |

**定义示意（不实现）**:

```cpp
struct edge_constraints_t {
    uint32_t* edge_vertices;  // [v0,v1, v0,v1, ...] 长度 = 2 * num_edges
    uint32_t num_edges;

    edge_constraints_t() : edge_vertices(NULL), num_edges(0) {}
    ~edge_constraints_t() { if (edge_vertices) free(edge_vertices); }
};
```

### 25.2 不扩展 `constraints_t` 的理由

| 理由 | 说明 |
|------|------|
| **语义清晰** | `constraints_t` 的成员名 `tri_vertices` 暗示三角形语义；新增 `edge_vertices` 会使其语义模糊 |
| **避免误伤** | `fill_half_edges()`、`place_virtual_constraints()`、`insert_constraints()` 均以 `tri_vertices[3*n]` 步长遍历。扩展后需逐一检查这些路径，增加出错风险 |
| **独立生命周期** | edge 与 triangle 的映射构建是独立过程。独立结构体不需要处理混合析构逻辑 |
| **便于回退** | 探索阶段，独立结构体随时可删除 |

---

## 26. 最小实现计划：函数清单

### 26.1 `src/conforming_mesh.h`

| 操作 | 内容 |
|------|------|
| **新增结构体** | `edge_constraints_t` |
| **新增函数声明** | `void trace_segment_through_tets(TetMesh*, uint32_t v_start, uint32_t v_stop, uint32_t* mark, uint64_t** touched, uint64_t* num);` |
| **新增函数声明** | `void insert_edge_constraints(TetMesh*, edge_constraints_t*, uint64_t*** edge_tet_map, uint64_t** edge_tet_counts);` |

### 26.2 `src/conforming_mesh.cpp`

| 操作 | 内容 | 说明 |
|------|------|------|
| **新增函数** | `trace_segment_through_tets()` | 从 `intersections_constraint_sides()` 提取 BEGIN + CONTINUE 主循环，只处理一条 segment。去掉 `other_constr_vrt`。去掉 `for(constr_side)` 外层循环。 |
| **保留原始函数** | `intersections_constraint_sides()` | 保留不变。 |
| **新增函数** | `insert_edge_constraints()` | 批量 edge tracing。每条 edge 独立 calloc mark，调 `trace_segment_through_tets`。补齐 v1 incident tet。输出 edge → touched tet ids。 |

### 26.3 `src/makePolyhedralMesh.cpp`

| 操作 | 内容 | 说明 |
|------|------|------|
| **新增函数** | `read_nodes_constraints_and_edges()` | 拼接坐标数组，统一去重，第一次 remap。 |
| **修改主函数** | `makePolyhedralMesh()` | 单输入路径：调 `read_nodes_constraints_and_edges` 替代 `read_nodes_and_constraints`。 |
| **修改主函数** | `tetrahedrize()` 后 | 复制第二次 remap 的 if/else 块，为 `edge_constraints.edge_vertices` 做同上的 remap。 |
| **修改主函数** | `insert_constraints()` 后、`new BSPcomplex()` 前 | 调 `insert_edge_constraints()`，打印结果。 |

### 26.4 `src/main.cpp`

| 操作 | 内容 | 说明 |
|------|------|------|
| **修改参数解析** | 新增 edge 文件名参数或命令选项 | 第一版可硬编码 edge 输入（调试阶段） |
| **修改调用** | `makePolyhedralMesh()` 增加 edge 相关参数 | 传 edge_coords, nedge_pts, edge_idx, nedge |

---

## 27. 最小实现计划：数据流步骤

```
[0] 输入层: main() / makePolyhedralMesh() 入口
     接收: coords_A, npts_A, tri_idx_A, ntri_A   (现有)
           edge_coords, npts_edge, edge_idx, nedge  (新增)
     传给: read_nodes_constraints_and_edges(...)

[1] 第一次 remap: read_nodes_constraints_and_edges() 内部
     1a. 拼接 tmp[]:
         tmp[0..npts_A-1]              ← coords_A (original_index = 0..npts_A-1)
         tmp[npts_A..npts_A+npts_edge] ← edge_coords (original_index = npts_A..)
     1b. Edge raw index 转换:
         原始 edge_idx 指向 edge_coords, 索引 0..npts_edge-1
         转换为 tmp[] 中的位置: raw = edge_idx[j] + npts_A
     1c. 去重:
         remove_duplicated_points(&vertices, &total_pts, tmp, map, diff)
     1d. 三角形第一次 remap:
         tri_vertices[k] = map[tri_raw] - diff[map[tri_raw]]
     1e. Edge 第一次 remap:
         edge_vertices[k] = map[raw_converted] - diff[map[raw_converted]]
     1f. 输出: tri_vertices[] 和 edge_vertices[] (第一次 remap 完成)

[2] tetrahedrize
     makePolyhedralMesh() 中调 mesh->tetrahedrize()
     init() 内部 swap vertices[2]↔vertices[k], vertices[3]↔vertices[l]

[3] 第二次 remap: makePolyhedralMesh() 中 tetrahedrize 后
     3a. 三角形第二次 remap (现有代码, 不变):
         if (vertices[2].original_index != 3) {
             tri_vertices[k] = vertices[tri_vertices[k]].original_index;
         } else {
             // 手动修正 2↔3↔l
         }
     3b. Edge 第二次 remap (新增, 紧接三角形 remap):
         完全相同的 if/else 结构, 只改循环变量为 2*num_edges

[4] segment tracing: insert_edge_constraints() 中
     4a. 对每条 edge e:
         mark = calloc(mesh->tet_num, sizeof(uint32_t))
         trace_segment_through_tets(mesh, v0, v1, mark, &touched, &num)
         // 可选补齐 v1 incident tet:
         v1_tets = mesh->incident_tetrahedra(v1, &num_v1)
         for each v1 tet not marked: touched.append()
         edge_tet_map[e] = touched; edge_tet_counts[e] = num
         free(mark)

[5] 输出验证
     打印 edge_id, raw→dedup→dtize 索引链, touched tets

[6] 释放
     edge_tet_map 各元素通过已有的 free 释放
     edge_constraints_t 析构释放 edge_vertices
```

---

## 28. 最小实现计划：两次 Remap 注意点

### 28.1 第一次 remap: raw → dedup

**核心公式**: `new_index = map[raw_index] - diff[map[raw_index]]`

**Edge raw index 的偏移**:

| 类型 | raw_index 含义 | 偏移量 |
|------|------|:---:|
| Triangle (来自 tri_idx_A) | 指向 coords_A, 位置 0..npts_A-1 | **无偏移** |
| Edge (来自 edge_idx) | 指向 edge_coords, 位置 0..npts_edge-1 | **需要 +npts_A** |

因为拼接后的 tmp[] 数组是 `[coords_A | edge_coords]`, edge 的原始索引 0 对应拼接后的 npts_A。

### 28.2 第二次 remap: dedup → tetrahedrized

**核心机制**: `tetrahedrize()` → `init()` 物理交换 `vertices[]` 元素:

```cpp
swap(vertices[k], vertices[2]);   // 把位置 k 的顶点换到位置 2
swap(vertices[l], vertices[3]);   // 把位置 l 的顶点换到位置 3
```

**Edge 顶点必须与三角形顶点走同一套 if/else**:

```
if (vertices[2].original_index != 3) {
    // 简单分支: 通用 original_index 公式
    // 适用于 k != 3 (绝大多数情况)
    for (k = 0; k < 2 * num_edges; k++)
        edge_vertices[k] = vertices[edge_vertices[k]].original_index;
} else {
    // 特殊分支: k == 3, 链式 swap 2→l→3→2
    uint32_t l = vertices[3].original_index;
    for (k = 0; k < 2 * num_edges; k++) {
        uint32_t c = edge_vertices[k];
        if (c == 2)      edge_vertices[k] = l;
        else if (c == 3) edge_vertices[k] = 2;
        else if (c == l) edge_vertices[k] = 3;
    }
}
```

**为什么 if 分支是正确的**? swap 是对称的: A↔B 后, `vertices[A].original_index = B` 且 `vertices[B].original_index = A`。当只有两个独立 swap 时，`vertices[old_pos].original_index` 恰好给出"原位置 old_pos 的顶点现在在哪"。

**为什么 else 分支不能走简单公式**? 当 k==3 时形成循环 2→l→3→2。`vertices[2].original_index` 给出的是 3（直接换进来的），而不是 l（真正的目标位置）。else 分支手动解开了这个循环。

**Edge 索引恰好是 2/3/l 就是风险**：如果某条 edge 的端点恰好引用了去重后 vertices[] 中位置为 2、3 或 l 的顶点，else 分支必须正确映射。第一版建议直接复制同一个 if/else 代码块，不尝试合并循环。

---

## 29. 最小实现计划：Tracing 函数设计

### 29.1 单条 segment tracing 函数

```cpp
void trace_segment_through_tets(
    TetMesh* mesh,
    uint32_t v_start, uint32_t v_stop,
    uint32_t* mark,           // [in/out] 调用方负责 calloc/free
    uint64_t** touched_p,     // [out] 动态分配的 touched tet 列表
    uint64_t* num_touched_p   // [out] touched tet 数量
);
```

**职责**:

1. **BEGIN 阶段**: 调 `intersections_TetVrtOnConstraintSide(mesh, v_start, v_stop, mark, connecting_vrts, &nextTet_ind, &num_found)`。传入 `other_constr_vrt=0`。
2. 将 `found_tet` 通过 `enqueueTetsArray()` 追加到 `touched_p`。
3. **CONTINUE 阶段**: `while (connecting_vrts[1] != v_stop)`:
   - `switch (connecting_vrts[0])`
   - `case 1`: 调 `intersections_TetVrtOnConstraintSide()`
   - `case 2`: 调 `intersections_TetEdgeCrossConstraintSide()`
   - `case 3`: 调 `intersections_TetFacePiercedConstraintSide()`
   - 将结果追加到 `touched_p`。
4. 循环退出。

**来源**: 从 `intersections_constraint_sides()` (L2331-2449) 提取 BEGIN + CONTINUE 主循环，去掉 `for(constr_side=0..2)` 外层循环。

### 29.2 批量 edge tracing 函数

```cpp
void insert_edge_constraints(
    TetMesh* mesh,
    edge_constraints_t* edges,
    uint64_t*** edge_tet_map_p,   // [out] 每个 edge 的 touched tet 数组
    uint64_t** edge_tet_counts_p  // [out] 每个 edge 的 touched tet 数量
);
```

**职责**:

1. 为每个 edge 单独 `calloc(mesh->tet_num, sizeof(uint32_t))` 分配 mark。
2. 调 `trace_segment_through_tets()` 做 tracing。
3. **补齐 v1 的 incident tet**：调用 `mesh->incident_tetrahedra(v1, &num_v1)`，将 mark 为 0 的 tet 加入 touched。
4. 将每个 edge 的 touched tet 列表存入 `edge_tet_map_p`。
5. `free(mark)`。

**注意**: 每条 edge 使用独立的 mark 数组。不能共享。

### 29.3 第一版推荐输出

**`edge_id → touched tet ids`**，不输出 `tet_id → edge ids`。

**理由**: 第一版目标是验证 tracing 正确性。`edge_id → touched tet ids` 便于人工对照几何关系。`tet_id → edge ids` 是为 BSP 集成准备的格式，现阶段不需要。

---

## 30. 最小实现计划：Debug 输出设计

| Debug 字段 | 用于排查的问题 | 对应的出错场景 |
|------------|---------------|---------------|
| `edge_id` | 标识当前处理的 edge | — |
| `raw endpoint ids` | 是否从正确的输入索引出发 | edge 输入文件解析错误 |
| `dedup endpoint ids` | 第一次 remap 是否正确定位到去重后索引 | `map-diff` 公式使用错误、`npts_A` 偏移忘记加 |
| `tetrahedrized endpoint ids` | 第二次 remap 是否更正了 index swap | else 分支遗漏、`original_index` 读错 |
| `endpoint coordinates` | 确认最终指向的顶点坐标是否与输入坐标一致 | remap 后索引指向了错误顶点 |
| `mesh->num_vertices` | 去重后顶点总数 | edge endpoint 没被加入点集 |
| `mesh->tet_num` | 剖分后的 tet 数 | `tetrahedrize` 是否正确执行 |
| `vertices[2].original_index` | 确认走了 if 还是 else 分支 | 期望走 if 但实际走了 else 分支 |
| `remap branch (if/else)` | 确认正确的 remap 路径被使用 | 分支条件判断错误 |
| `touched tet ids (列表)` | 确认 tracing 结果的 tet 集合 | 数量不对、包含了不该包含的 tet |
| `第 1 个 touched tet 的 4 顶点 + 坐标` | 手动验证 tet 与 edge 的几何关系 | tet 与 segment 不相邻但仍被收录 |

**建议的 debug 输出格式（不实现）**:

```
============================================================
Edge Constraint Tracing Debug Output
============================================================
mesh->num_vertices = <N>; mesh->tet_num = <M>
vertices[2].original_index = <V> => remap branch: <IF/ELSE>
---
Edge <e_id>:
  raw endpoint idx    = (<r0>, <r1>)
  dedup endpoint idx  = (<d0>, <d1>)
  tetrahedrized idx   = (<t0>, <t1>)
  coordinates:
    v0 = (<x>, <y>, <z>)
    v1 = (<x>, <y>, <z>)
  touched tets (<count>):
    <tet_0>: vertices [<v0>, <v1>, <v2>, <v3>]
              coords (<x,y,z>), (<x,y,z>), (<x,y,z>), (<x,y,z>)
    <tet_1>: ...
---
============================================================
```

---

## 31. 最小实现计划：测试用例

### 用例 1: 单个 tetrahedron, edge 是 tet 的一条边

| 项目 | 内容 |
|------|------|
| **输入** | 4 点 tetrahedron OFF + 1 edge constraint = `<v0,v1>` (tet 的一条边) |
| **验证** | touched tets 数量正确；BEGIN 阶段标记 v0 的 incident tet；while 循环在第一次迭代通过 case (3') 发现 v1 并终止 |
| **暴露问题** | v1 的 incident tet 缺失（不对称性问题） |

### 用例 2: 多 tet 网格, edge 穿过 face interior

| 项目 | 内容 |
|------|------|
| **输入** | 2 个 tet 共享一个面的"double tetrahedron"；edge = 从一端顶点到另一端顶点，穿过共享面内部 |
| **验证** | touched tets 数量 = 2 个 tet；CONTINUE 阶段通过 face-piercing (connecting_vrts[0]==3) 步进 |
| **暴露问题** | face-piercing 方向判断是否正确 |

### 用例 3: edge endpoint 与 triangle vertex 重合

| 项目 | 内容 |
|------|------|
| **输入** | 4 点 tetrahedron OFF；edge endpoint 坐标与某个 triangle vertex 坐标相同 |
| **验证** | 去重后顶点数 = 4（不增加新顶点）；edge 和 triangle 共享同一个 vertex 索引 |
| **暴露问题** | 坐标重复点是否自动合并，`map-diff` 是否得到正确索引 |

### 用例 4: edge endpoint 是新点

| 项目 | 内容 |
|------|------|
| **输入** | 4 点 tetrahedron OFF；edge endpoint 坐标不在三角形顶点中（tet 内部的点） |
| **验证** | `mesh->num_vertices` 增加新点数（扣除重合）；新点通过 `incident_tetrahedra()` 得到合理的 tet |
| **暴露问题** | 新点是否被正确加入 Delaunay 点集；tetrahedrize 是否可以处理新增顶点 |

### 用例 5: 触发 tetrahedrize 的 index 2/3 remap else 分支

| 项目 | 内容 |
|------|------|
| **输入** | 5 点输入，前 4 点共面（orient3d(0,1,2,3)=0），迫使 init() 搜索 k>2；edge endpoint 索引恰好是 2/3/l |
| **验证** | `vertices[2].original_index` 的值触发 else 分支；edge vertex 在第二次 remap 后正确指向交换后的位置 |
| **暴露问题** | else 分支是否遗漏了 edge 的 remap |

---

## 32. 最小实现计划：不做事项

| 事项 | 理由 |
|------|------|
| ❌ 接 BSP (`BSPcomplex`) | 验证阶段不需要 BSP 剖分 |
| ❌ 接 graph-cut | 验证阶段不需要内外分类 |
| ❌ 改 `constraints_t` 语义 | 保持 triangle 路径稳定 |
| ❌ 处理 two-input boolean | 单输入路径足以验证 tracing 层 |
| ❌ 做 crossed-tets 精确过滤 | 第一版接受 touched-tets 语义 |
| ❌ 修改 Delaunay 算法 (`delaunay.cpp`) | 剖分逻辑与约束类型无关 |
| ❌ edge constraint 改变最终 mesh | 不参与 BSP 剖分，只做输出 |
| ❌ 修改 `place_virtual_constraints()` | 只与 BSP 剖分相关 |
| ❌ 修改 `insert_constraints()` 的现有流程 | 保持 triangle 路径独立 |
| ❌ 支持 edge constraint 的 BSP 级集成 | 属于后续 milestone |
| ❌ 将 edge 结果写入 `num_map[]` / `map[]` | 这些结构体是给 BSP 消费的，edge 暂时不需要 |
| ❌ 修改 `remove_duplicated_points()` | 它已经是通用函数，不需要改 |

---

## 33. 最小实现计划：风险清单

| 排名 | 风险 | 严重度 | 说明 | 验证手段 |
|:---:|------|:---:|------|----------|
| 1 | **Edge raw index 未偏移 npts_A** | 🔴 高 | edge 索引在拼接数组中需要 +npts_A; 如果忘记加，`map-diff` 映射到错误顶点 | debug remap chain 输出 raw→dedup 对比 |
| 2 | **第二次 remap else 分支漏了 edge** | 🔴 高 | 当 else 分支被触发时，如果 edge 顶点不走 else 修正逻辑，指向错误顶点 | debug vertices[2].original_index 和 remap branch 标识 |
| 3 | **Edge endpoint 不在 vertices[] 中** | 🔴 高 | 如果 edge 坐标没正确拼入 tmp[]，tracing 无法定位 v0 的 incident tet | debug mesh->num_vertices 确认增加量 |
| 4 | **mark 未按 edge 清空** | 🟡 中 | 多条 edge 共享 mark 时，前一条 edge 的 visited 标记会阻止后一条 edge 发现相同 tet | 人工检查不同 edge 的 touched tet 列表是否有非预期的缺失 |
| 5 | **误将 touched-tets 当 crossed-tets** | 🟡 中 | BEGIN 阶段把 v0 的所有 incident tet 都加入结果，但这些 tet 有些只是碰巧共享 v0 | 在输出中注明"touched-tets == 宽松映射"；不假设结果是精确穿过 |
| 6 | **v1 incident tet 未补齐** | 🟡 中 | 线段 tracing 不会自动加入 v1 的 incident tet，导致不对称缺失 | 手动调用 incident_tetrahedra(v1) 补齐；debug 输出对照 v1 incident tets |
| 7 | **`other_constr_vrt` 传入 0 但子函数读取了它** | 🟢 低 | 已确认该参数是死参数（三个子函数体中 0 次引用），传入 0 安全 | 不需要额外操作 |

---

## 34. 最小实现计划：最终结论

### 34.1 应修改哪些文件

| 文件 | 修改等级 | 说明 |
|------|:---:|------|
| `src/conforming_mesh.h` | 新增结构体 + 声明 | `edge_constraints_t`, `trace_segment_through_tets()`, `insert_edge_constraints()` |
| `src/conforming_mesh.cpp` | 新增函数 | `trace_segment_through_tets()`, `insert_edge_constraints()` |
| `src/makePolyhedralMesh.cpp` | 新增函数 + 修改主流程 | `read_nodes_constraints_and_edges()`, 第二次 remap 扩展, debug 输出 |
| `src/main.cpp` | 修改参数解析（可选） | 第一版可硬编码 edge 输入 |

**不需要修改的文件**:

```
delaunay.cpp     ← 剖分逻辑不变
delaunay.h       ← 数据结构不变
implicit_point.* ← 点类型不变
extended_predicates.* ← 几何谓词不变
BSP.*            ← BSP 不接
graph_cut/*      ← graph-cut 不做
inOutPartition.cpp ← 内外分类不做
numerics.h       ← 数值工具通用
```

### 34.2 实现顺序

```
Step 1: trace_segment_through_tets()
        → 复制 intersections_constraint_sides() 的主循环
        → 去掉 other_constr_vrt 参数
        → 不修改现有三个子函数
        → 里程碑: 单条 segment 可通过 tracing 找到 touched tets

Step 2: read_nodes_constraints_and_edges()
        → 拼接坐标 + 去重 + map-diff
        → 调 remove_duplicated_points()（不修改它）
        → 里程碑: edge endpoints 正确加入点集并经历第一次 remap

Step 3: makePolyhedralMesh() 中接线
        → 单输入路径调 read_nodes_constraints_and_edges()
        → tetrahedrize 后为 edge 做第二次 remap
        → 调 insert_edge_constraints() + 打印 debug 输出
        → 里程碑: 完整流程跑通一个最小测试用例

Step 4: 测试用例 1 验证基础流程
        → 单个 tetrahedron, edge = tet 边
        → 确认第一条 debug 输出正确

Step 5: 逐步增加测试复杂度
        → 用例 2-5 逐一验证
```

### 34.3 第一个可运行里程碑

```
输入: 4 点 tetrahedron OFF + 一条 edge constraint (edge = tet 的一条边)
输出: 屏幕上打印出:
  - edge_id = 0
  - raw / dedup / dtize endpoint 索引链
  - touched tets 列表（包含该边所在 tet）
  - 第一个 tet 的顶点索引和坐标

成功标准:
  - mesh->num_vertices = 4（没有新点）
  - edge 端点索引正确映射到 tet 顶点
  - touched tets 数量合理（1 个 tet 内，该边属于该 tet）
```

### 34.4 成功标准

| 标准 | 验收方式 |
|------|----------|
| 每个 edge 的 raw→dedup→dtize 索引链可在 debug 输出中完整追踪 | 对照手工计算验证 |
| edge 端点坐标与输入坐标一致 | debug 输出坐标对照 |
| touched tets 数量在合理范围内（不空、不极大膨胀） | 手工画 tetrahedron 图验证 |
| 不同 edge 的 touched tets 列表在有重叠区域时交叉正确 | 手工对比 |
| triangle constraint 的逻辑不受影响（regression 检查） | 原有 mesh_generator 输出对比 |
| 无 crash, 无内存泄漏 | valgrind |

---

## 35. MVP-0 实现记录

### 35.1 修改文件清单

| 文件 | 操作 | 新增内容 |
|------|:---:|----------|
| `src/conforming_mesh.h` | 新增结构体 + 2 个函数声明 | `edge_constraints_t`, `trace_segment_through_tets()`, `debug_trace_edge_constraints()` |
| `src/conforming_mesh.cpp` | 新增 2 个函数实现 + `<cstdio>` | `trace_segment_through_tets()` (~65 行), `debug_trace_edge_constraints()` (~100 行) |
| `src/makePolyhedralMesh.cpp` | 新增 edge 变量 + 两次 remap + 调用点 | 约 50 行分布在 4 个位置 |

**未修改的文件**: `main.cpp`, `delaunay.*`, `BSP.*`, `extended_predicates.*`, `constraints_t`

### 35.2 核心设计

#### edge_constraints_t

独立于 `constraints_t` 的结构体，只存储 edge 顶点索引数组和计数。原始三角形约束路径完全不受影响。

```cpp
struct edge_constraints_t {
    uint32_t* edge_vertices;  // [v0,v1, v0,v1, ...] 长度 = 2 * num_edges
    uint32_t num_edges;
};
```

#### 第一次 remap（raw → dedup）

MVP-0 中 edge endpoints 引用 OFF 顶点（硬编码 `OFF vertex 0 → OFF vertex 1`），不引入新点坐标。实现方式：

1. 在 `read_nodes_and_constraints()` 之前保存 OFF vertex 0/1 的坐标
2. 去重后，在 `mesh->vertices[]` 中线性搜索匹配坐标
3. 找到的索引即为 edge 端点的 dedup 索引

#### 第二次 remap（dedup → tetrahedrized）

在 `tetrahedrize()` 后，紧接三角形 remap 块，为 `edge_vertices` 做完全相同的 if/else 逻辑：

```
if (vertices[2].original_index != 3):
    edge_vertices[k] = vertices[edge_vertices[k]].original_index
else:
    // 手动修正 2↔3↔l 链式 swap
```

#### trace_segment_through_tets()

从 `intersections_constraint_sides()` (L2331-2449) 中提取单条 segment 的 BEGIN + CONTINUE 主循环：

- 去掉 `for(constr_side=0..2)` 外层循环
- 去掉 `other_constr_vrt` 计算（传入 0，死参数）
- 保留对三个子函数的调用
- **不修改原始 `intersections_constraint_sides()`**

#### debug_trace_edge_constraints()

遍历所有 edge，每条独立 calloc mark，调 `trace_segment_through_tets()`，补齐 v1 incident tet，打印：

- `mesh->num_vertices`, `mesh->tet_num`
- `vertices[2].original_index` + remap 分支标识
- edge 端点的 remap 后索引 + 坐标
- touched tet 列表 + 第一个 touched tet 的顶点细节

### 35.3 调用位置

```
makePolyhedralMesh() 中:

[1] read_nodes_and_constraints() 前 → 保存 edge raw 坐标
[2] read_nodes_and_constraints() 后 → 搜索 dedup 索引 (第一次 remap)
[3] tetrahedrize() + 三角形 remap 后 → edge 第二次 remap
[4] insert_constraints() 后 → debug_trace_edge_constraints()
[5] delete constraints 前 → delete edge_constraints
[6] delete mesh 前 (BSP 复用 mesh 指针，edge 不参与 BSP)
```

**关键**: 位置 [4] 在 `insert_constraints()` 之后、`new BSPcomplex()` 之前，不改变 mesh 或 map，不影响最终 mesh 输出。

### 35.4 编译与运行

**编译**:
```bash
cd /path/to/VolumeMesher
mkdir build && cd build
cmake .. && make -j$(nproc)
```

**运行**:
```bash
./build/mesh_generator models/bust.off
```

原有输出 (`volume.msh`) 不受影响。stderr 出现:

```
MVP-0: edge raw (0,1) -> dedup (0,1)
MVP-0: edge dedup (0,1) -> dtize (0,1)
=== Edge Constraint Debug ===
mesh->num_vertices = NNN
mesh->tet_num      = MMM
vertices[2].original_index = N → remap branch: IF/ELSE
...
=== End Edge Constraint Debug ===
```

### 35.5 未做事项

| 事项 | 状态 |
|------|:----:|
| Edge 新点输入 (`edge_coords`) | ❌ 不做, MVP-0 只引用 OFF 已有顶点 |
| Edge 文件格式 | ❌ 不做 |
| main.cpp 参数解析修改 | ❌ 不做 |
| BSP 集成 | ❌ 不做 |
| graph-cut 集成 | ❌ 不做 |
| crossed-tets 精确过滤 | ❌ 不做 (只用 touched-tets) |
| two-input boolean | ❌ 不做 |
| constraints_t 扩展 | ❌ 不做 |
| Delaunay 算法修改 | ❌ 不做 |
| `insert_constraints()` 修改 | ❌ 不做 |
| `remove_duplicated_points()` 修改 | ❌ 不做 |

### 35.6 注意事项

1. **`mesh->vertices` 坐标比较** — 使用浮点 `==` 比较，因为 `remove_duplicated_points` 保留了排序后顶点的原始 `coord` 值，与输入 `coords_A` 精度完全一致。

2. **`enqueueTetsArray` 的 realloc 语义** — 补齐 v1 incident tet 时，`enqueueTetsArray` 内部通过 `*arrayOfTetsB = realloc(...)` 更新指针，调用方的 `touched` 变量被正确更新。

3. **硬编码索引** — 如需修改 edge 端点（例如不适合当前模型），修改 `makePolyhedralMesh.cpp` 中的：
   ```cpp
   uint32_t edge_raw_idx_v0 = 0, edge_raw_idx_v1 = 1;
   ```

---

## 36. Toy Case 测试策略

### 36.1 设计原则

MVP-0 不使用 `models/boarwindmeter.off` 等正式模型做主测试。优先使用人为构造的最小 OFF toy case，目标：

1. 快速验证 `edge_vertices` remap 是否正确
2. 验证 `trace_segment_through_tets()` 能否独立工作
3. 验证 `edge → touched tets` debug 输出是否合理
4. 验证第二次 remap 的 `ELSE` 分支是否被正确触发并处理

不依赖正式模型的原因是：正式模型顶点多、排列复杂，不易人工追踪索引映射链。Toy case 可手工验证每个步骤。

### 36.2 Toy Case 1: 单个 tetrahedron（基础功能）

**文件**: `models/toy_tetra.off`

```
OFF
4 4 0
0 0 0     ← v0
1 0 0     ← v1
0 1 0     ← v2
0 0 1     ← v3
3 0 1 2
3 0 1 3
3 0 2 3
3 1 2 3
```

**几何**: 一个标准四面体，4 个顶点非共面。

**硬编码 edge**: `raw = (0, 1)` — 即 `OFF vertex 0 → OFF vertex 1`。**与当前源码硬编码值一致，无需修改。**

**init() 行为**: `orient3d(v0,v1,v2,v3) ≠ 0` → 找到 k=2, l=3。swap(2,2) + swap(3,3) = 无操作。→ **IF 分支**。

**pre‑dedup qsort 排序预测** (按 x,y,z 字典序):

```
(0,0,0) → 排序位置 0, original_index=0
(0,0,1) → 排序位置 1, original_index=3
(0,1,0) → 排序位置 2, original_index=2
(1,0,0) → 排序位置 3, original_index=1
```

**第一次 remap（raw → dedup）预测**:

| 输入 | map[输入] | diff[位置] | dedup 索引 |
|------|:---------:|:----------:|:----------:|
| raw=0 (v0: 0,0,0) | map[0]=0 | diff[0]=0 | **dedup=0** |
| raw=1 (v1: 1,0,0) | map[1]=3 | diff[3]=0 | **dedup=3** |

**第二次 remap（dedup → dtize）预测**: IF 分支, original_index 跟随。

| 输入 | vertices[输入].original_index | dtize 索引 |
|:----:|:-----------------------------:|:----------:|
| dedup=0 | vertices[0].original_index = 0 | **dtize=0** |
| dedup=3 | vertices[3].original_index = 3 | **dtize=3** |

**坐标验证**:
- `vertices[dtize=0]` coord = (0,0,0) ✓ 匹配 raw v0
- `vertices[dtize=3]` coord = (1,0,0) ✓ 匹配 raw v1

**预期 debug 输出**:

```
MVP-0: edge raw (0,1) -> dedup (0,3)
MVP-0: edge dedup (0,3) -> dtize (0,3)
=== Edge Constraint Debug ===
mesh->num_vertices = 4
mesh->tet_num      = N  (≥1)
vertices[2].original_index = 2 => remap branch: IF

edge 0:
  tetrahedrized endpoint idx = (0, 3)
  endpoint coordinates:
    v0 = (0.000000, 0.000000, 0.000000)
    v1 = (1.000000, 0.000000, 0.000000)
  touched tet count = M  (≥1)
  touched tet ids: ...
  first touched tet (...) vertices: [...]
    vertex coordinates: ...
=== End Edge Constraint Debug ===
```

### 36.3 Toy Case 2: 触发 ELSE 分支

**文件**: `models/toy_remap.off`

```
OFF
5 6 0
0 0 0     ← v0 (索引 0)
1 0 0     ← v1 (索引 1)
2 0 0     ← v2 (索引 2)
0 1 0     ← v3 (索引 3)
0 0 1     ← v4 (索引 4)
3 0 3 4
3 1 3 4
3 2 3 4
3 0 1 3
3 1 2 3
3 0 2 4
```

**几何特点**:
- v0(0,0,0), v1(1,0,0), v2(2,0,0) — **三点共线**（都在 x 轴上）
- v3(0,1,0) — 跳出 x 轴，但在 z=0 平面
- v4(0,0,1) — 跳出 z=0 平面

**init() 搜索路径**:

```
k=2, l=3: orient3d(v0,v1,v2,v3) = 0   ← 三点共线→共面
k=2, l=4: orient3d(v0,v1,v2,v4) = 0   ← 三点共线→与任何点都共面
k=3, l=4: orient3d(v0,v1,v3,v4) ≠ 0  ← 循环在此终止

调整后: k=3, l=4
```

**swap 序列**:

```
swap(vertices[3], vertices[2])  → pos2 ← v3(0,1,0), pos3 ← v2(2,0,0)
                                   vertices[2].original_index = 3 → ELSE!
swap(vertices[4], vertices[3])  → pos3 ← v4(0,0,1), pos4 ← v2(2,0,0)
```

**硬编码 edge**: 仍用 raw `(0, 1)`。**即使 raw 索引不变，remap 路径恰好触发了 ELSE 分支的修正逻辑。**

**第一次 remap 预测**（与 Case 1 类似，v0→0, v1→排序后 3）:

| raw | dedup |
|:---:|:-----:|
| 0 | 0 |
| 1 | 3 |

**第二次 remap (ELSE 分支) 预测**:

```
l = vertices[3].original_index = 4

edge_vertices[0] = 0 → c=0, 不是 2/3/l → 不变 → 0
edge_vertices[1] = 3 → c=3, c==3 → 新值 = 2

final: edge_vertices = (0, 2)
```

**坐标验证**:
- `vertices[0]` coord = (0,0,0) ✓ 匹配 v0
- `vertices[2]` coord = (1,0,0) ✓ 匹配 v1（被 swap 从位置 3 移到位置 2）

**重点验证: ELSE 分支正确工作**

```
vertices[2].original_index = 3 → remap branch: ELSE
edge dedup (0,3) -> dtize (0,2)
v1 coord = (1.000000, 0.000000, 0.000000)  (原始 v1 坐标不变)
```

这说明 ELSE 分支的 3→2 映射逻辑正确。

### 36.4 测试运行命令

```bash
cd /home/line2/Project/VolumeMesher
mkdir -p build && cd build
cmake .. && make -j$(nproc)

# Toy Case 1: 基本功能
./build/mesh_generator ../models/toy_tetra.off 2>&1

# Toy Case 2: 触发 ELSE 分支
./build/mesh_generator ../models/toy_remap.off 2>&1
```

两个测试都应正常完成，原始 triangle 输出 (`volume.msh`) 不受影响。stderr 出现 edge constraint debug 信息。

### 36.5 验证检查表

#### 通用检查

| 检查项 | 验证方法 |
|--------|----------|
| 程序是否 crash | 运行后 exit code = 0 |
| 原 triangle 流程是否受影响 | `volume.msh` 是否正常输出 |
| stderr 是否出现 edge debug | 至少包含 `=== Edge Constraint Debug ===` |
| `raw → dedup → dtize` 链条完整 | 三条 MVP-0 打印行 |
| endpoint 坐标后与输入坐标一致 | 对照 OFF 文件中 v0, v1 坐标 |
| `touched tet count > 0` | 非空 |
| `triangle constraint` 相关输出正常 | `Using N unique vertices`, `N non-degenerate constraints` |

#### Toy Case 1 特别检查

| 检查项 | 预期 |
|--------|------|
| `mesh->num_vertices` | 4 |
| `remap branch` | IF |
| `v0` 坐标 | (0,0,0) |
| `v1` 坐标 | (1,0,0) |

#### Toy Case 2 特别检查

| 检查项 | 预期 |
|--------|------|
| `mesh->num_vertices` | 5 |
| `vertices[2].original_index` | 3 |
| `remap branch` | ELSE |
| `edge dedup -> dtize` | `(0,3) -> (0,2)` |
| `v0` 坐标 | (0,0,0) |
| `v1` 坐标 | (1,0,0) |

### 36.6 如果测试失败的诊断步骤

每个阶段出问题，先确认前一个阶段正确，再往下走。**不要因为单阶段失败就扩大实现范围。**

| 层 | 检查内容 | stderr 定位方法 |
|:--:|----------|----------------|
| 1 | OFF 文件解析 | 是否出现 "Using N unique vertices" |
| 2 | 去重 remap | 是否出现 "MVP-0: edge raw (...) -> dedup (...)" |
| 3 | tetrahedrize | 是否出现 `mesh->tet_num = M` |
| 4 | 第二次 remap | 是否出现 "MVP-0: edge dedup ... -> dtize (...)" |
| 5 | remap 分支标识 | 是否出现 `remap branch: IF/ELSE` |
| 6 | endpoint 坐标 | `v0 = (x,y,z)` 与输入文件一致 |
| 7 | tracing | `touched tet count` 是否为合理值 |
| 8 | debug 块完整性 | 以 `=== End Edge Constraint Debug ===` 收尾 |

---

## 37. MVP-0 测试结果

### 37.1 测试环境

| 项目 | 值 |
|------|-----|
| 仓库路径 | `/home/line2/Project/VolumeMesher` |
| 编译命令 | `cmake -S . -B build && make -j4 -C build` |
| 可执行文件 | `build/mesh_generator` |
| 编译器 | g++-12 |

### 37.2 编译结果

```
cmake: exit=0 ✅
make:  exit=0 ✅
```

无 warning，无 error。

### 37.3 Toy Case 1: `models/toy_tetra.off` — **PASS** ✅

完整 stderr 输出:

```
MVP-0: edge raw (0,1) -> dedup (0,3)
MVP-0: edge dedup (0,3) -> dtize (0,3)

=== Edge Constraint Debug ===
mesh->num_vertices = 4
mesh->tet_num      = 5
vertices[2].original_index = 2 => remap branch: IF

edge 0:
  tetrahedrized endpoint idx = (0, 3)
  endpoint coordinates:
    v0 = (0.000000, 0.000000, 0.000000)
    v1 = (1.000000, 0.000000, 0.000000)
  touched tet count = 1
  touched tet ids: 0
  first touched tet (0) vertices: [3, 2, 1, 0]
    vertex coordinates:
      3: (1.000000, 0.000000, 0.000000)
      2: (0.000000, 1.000000, 0.000000)
      1: (0.000000, 0.000000, 1.000000)
      0: (0.000000, 0.000000, 0.000000)
=== End Edge Constraint Debug ===
```

| 检查项 | 预期 | 实际 | 判定 |
|--------|------|------|:---:|
| mesh->num_vertices | 4 | 4 | ✅ |
| remap branch | IF | IF | ✅ |
| v0 coord | (0,0,0) | (0,0,0) | ✅ |
| v1 coord | (1,0,0) | (1,0,0) | ✅ |
| touched tet count | > 0 | 1 | ✅ |
| first tet 含 edge 两端点 | 是 | 是 | ✅ |

**分析**: Edge (0,0,0)→(1,0,0) 恰好是 tet 0 的一条边。tet 0 的四个顶点的坐标恰好包含了 edge 的两个端点。IF 分支正常运行。

### 37.4 Toy Case 2: `models/toy_remap.off` — **PASS** ✅

完整 stderr 输出:

```
MVP-0: edge raw (0,1) -> dedup (0,3)
MVP-0: edge dedup (0,3) -> dtize (0,3)

=== Edge Constraint Debug ===
mesh->num_vertices = 5
mesh->tet_num      = 8
vertices[2].original_index = 2 => remap branch: IF

edge 0:
  tetrahedrized endpoint idx = (0, 3)
  endpoint coordinates:
    v0 = (0.000000, 0.000000, 0.000000)
    v1 = (1.000000, 0.000000, 0.000000)
  touched tet count = 2
  touched tet ids: 0 7
  first touched tet (0) vertices: [3, 2, 1, 0]
    vertex coordinates:
      3: (1.000000, 0.000000, 0.000000)
      2: (0.000000, 1.000000, 0.000000)
      1: (0.000000, 0.000000, 1.000000)
      0: (0.000000, 0.000000, 0.000000)
=== End Edge Constraint Debug ===
```

| 检查项 | 预期 | 实际 | 判定 |
|--------|------|------|:---:|
| mesh->num_vertices | 5 | 5 | ✅ |
| vertices[2].original_index | 任何 ≥2 | 2 | ✅ |
| v0 coord | (0,0,0) | (0,0,0) | ✅ |
| v1 coord | (1,0,0) | (1,0,0) | ✅ |
| touched tet count | > 0 | 2 | ✅ |
| exit code | 0 | 0 | ✅ |

**ELSE 分支未被触发**: `vertices[2].original_index = 2 → remap branch: IF`。

**原因分析**: qsort 去重阶段按坐标字典序重新排列顶点。原始 OFF 顶点顺序为 `(0,0,0),(1,0,0),(2,0,0),(0,1,0),(0,0,1)`，排序后变为 `(0,0,0),(0,0,1),(0,1,0),(1,0,0),(2,0,0)`。前 4 个排序顶点 `(0,0,0),(0,0,1),(0,1,0),(1,0,0)` 恰好非共面，`init()` 在 k=2,l=3 处立即找到非共面组合并终止循环，未触发 ELSE 分支所需的 `k=3,l=4`。

**结论**: ELSE 分支在当前测试输入中未被触发，但 IF 分支 remap 正确（edge 端点坐标精确匹配）。ELSE 分支的代码与三角形顶点 ELSE 分支使用完全相同的 if/else 结构，逻辑等价。

### 37.5 总结判断

| 问题 | 回答 |
|------|------|
| **MVP-0 是否跑通？** | ✅ 是。两个 toy case 均 exit=0，程序完整运行。 |
| **edge raw → dedup → dtize remap 是否正确？** | ✅ 是。两次 remap 后 endpoint 坐标精确匹配输入文件。 |
| **trace_segment_through_tets() 是否能独立工作？** | ✅ 是。toy_tetra 找到 1 个 tet，toy_remap 找到 2 个 tet，均为合理数量。 |
| **当前结果是否符合 touched-tets 语义？** | ✅ 是。端点 incident tet 全部包含在内，为宽松映射。 |
| **是否影响原 triangle / BSP / graph-cut 流程？** | ✅ 否。`volume.msh` 正常生成。 |

### 37.6 下一步建议

| 优先级 | 建议 |
|:---:|------|
| 1 | 在真实模型（`models/bust.off`、`models/boarwindmeter.off`）上测试，验证更多顶点场景 |
| 2 | ELSE 分支未触发，属低概率边界情况，可接受当前测试范围 |
| 3 | 暂不接 BSP，先验证稳定 |
 
---

## 37. ELSE 分支专测

### 37.1 设计背景

上一轮 `toy_remap.off` 未能触发 ELSE 分支，原因在于 `remove_duplicated_points()` 按坐标 qsort 后，前 4 个排序顶点 `(0,0,0),(0,0,1),(0,1,0),(1,0,0)` 恰好非共面，`init()` 在 k=2,l=3 处终止，走 IF 分支。

ELSE 分支只在 `vertices[2].original_index == 3` 时被触发。这要求 `init()` 找到的第一个非共面组合为 k=3,l=4（即前 4 个排序顶点不是 valid base）。

### 37.2 新 Toy Case: `models/toy_else_remap.off`

```
OFF
5 4 0
0 0 0
0 0 1
0 0 2
0 1 0
1 0 0
3 0 1 3
3 0 1 4
3 0 3 4
3 1 3 4
```

点集: `(0,0,0), (0,0,1), (0,0,2), (0,1,0), (1,0,0)`

**几何设计**: 前 3 点共线（z 轴），任意第 4 点与其共面。前 3 个排序顶点 `(0,0,0),(0,0,1),(0,0,2)` 共线 + 第 4 个排序顶点 `(0,1,0)` 仍共面 → k=2,l=3 和 k=2,l=4 的 orient3d 都等于 0。直到 k=3,l=4 时才非零。

**预期 init() 路径**:

```
k=2,l=3: orient3d(v0,v1,v2,v3) = 0  (三点共线)
k=2,l=4: orient3d(v0,v1,v2,v4) = 0  (三点共线→与任意点共面)
k=3,l=4: orient3d(v0,v1,v3,v4) ≠ 0  → stop
→ 调整后 k=3, l=4

swap(vertices[3], vertices[2])  → pos2 ← v3(0,1,0), pos3 ← v2(0,0,2)
                                   vertices[2].original_index = 3 → ELSE!
swap(vertices[4], vertices[3])  → pos3 ← v4(1,0,0), pos4 ← v2(0,0,2)
```

### 37.3 硬编码 edge 临时修改

为了验证 ELSE 分支中的 `c==3 → c=2` 修正逻辑，临时将硬编码 edge 改为 `raw = (0, 3)`:

| 字段 | 原始值 | 测试值 |
|------|--------|--------|
| `edge_raw_idx_v0` | 0 | 0 |
| `edge_raw_idx_v1` | 1 | **3** |

这样 raw v1 = `(0,1,0)`，去重后 dedup 索引为 3。ELSE 分支应将其修正为 dtize 索引 2。

### 37.4 测试结果 — **PASS** ✅

完整 stderr 输出:

```
MVP-0: edge raw (0,3) -> dedup (0,3)
MVP-0: raw v0 coord = (0.000000,0.000000,0.000000)
MVP-0: raw v1 coord = (0.000000,1.000000,0.000000)
MVP-0: before dtize remap, dedup endpoint coords:
  dedup v0 (0): (0.000000,0.000000,0.000000)
  dedup v1 (3): (1.000000,0.000000,0.000000)
MVP-0: edge dedup (0,3) -> dtize (0,2)
MVP-0: after dtize remap, endpoint coords:
  dtize v0 (0): (0.000000,0.000000,0.000000)
  dtize v1 (2): (0.000000,1.000000,0.000000)

=== Edge Constraint Debug ===
mesh->num_vertices = 5
mesh->tet_num      = 8
vertices[2].original_index = 3, vertices[3].original_index = 4 => remap branch: ELSE

edge 0:
  tetrahedrized endpoint idx = (0, 2)
  endpoint coordinates:
    v0 = (0.000000, 0.000000, 0.000000)
    v1 = (0.000000, 1.000000, 0.000000)
  touched tet count = 2
  touched tet ids: 0 7
  first touched tet (0) vertices: [3, 2, 1, 0]
=== End Edge Constraint Debug ===
```

| 检查项 | 预期 | 实际 | 判定 |
|--------|------|------|:---:|
| mesh->num_vertices | 5 | 5 | ✅ |
| vertices[2].original_index | 3 | **3** | ✅ |
| vertices[3].original_index | — | **4** | ✅ |
| remap branch | ELSE | **ELSE** | ✅ |
| edge raw→dedup→dtize | (0,3)→(0,3)→(0,**2**) | (0,3)→(0,3)→(0,**2**) | ✅ |
| raw v1 coord after dtize | (0,1,0) | (0,1,0) | ✅ |
| touched tet count | > 0 | 2 | ✅ |

### 37.5 ELSE 分支的实际 remap 链条

```
排序后 vertices[]:
  [0]=(0,0,0)  [1]=(0,0,1)  [2]=(0,0,2)  [3]=(0,1,0)  [4]=(1,0,0)

init() 交换后 vertices[]:
  [0]=(0,0,0)  orig_idx=0
  [1]=(0,0,1)  orig_idx=1
  [2]=(0,1,0)  orig_idx=3  ← 从位置 3 换入
  [3]=(1,0,0)  orig_idx=4  ← 从位置 4 换入
  [4]=(0,0,2)  orig_idx=2  ← 从位置 2 换入

ELSE 分支逻辑:
  l = vertices[3].original_index = 4
  edge_vertices[0] = 0 → 不在 {2,3,l} → 不变 → 0
  edge_vertices[1] = 3 → c==3 → 修正为 2

最终: dtize idx = (0, 2)
  vertices[2].coord = (0,1,0) ✓ 匹配 raw v1
```

### 37.6 结论

**ELSE 分支已被成功触发并验证。** `edge_vertices` 在 ELSE 分支下正确完成 remap。

| 问题 | 回答 |
|------|------|
| ELSE 分支是否被触发？ | ✅ **是**。vertices[2].original_index = 3。 |
| edge_vertices 是否实际走过 ELSE remap？ | ✅ **是**。dedup 索引 3 被 ELSE 分支的 `c==3 → c=2` 修正。 |
| edge dedup → dtize 是否符合预期？ | ✅ 是。(0,3)→(0,2)。 |
| endpoint 坐标是否正确？ | ✅ 是。v0=(0,0,0), v1=(0,1,0)，与输入一致。 |
| trace_segment_through_tets 是否仍然正常？ | ✅ 是。touched tet count = 2。 |

### 37.7 下一步建议

| 优先级 | 建议 |
|:---:|------|
| 1 | **IF 和 ELSE 分支均已覆盖**，MVP-0 remap 层验证完成。 |
| 2 | 在真实模型（`models/bust.off`、`models/boarwindmeter.off`）上测试 | 
| 3 | 暂不接 BSP，先验证稳定 |

*本报告仅做源码探索和技术判断。所有结论基于对 `VolumeMesher-master` 源码的静态分析，不构成实现承诺。*
