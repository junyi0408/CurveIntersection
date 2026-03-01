# CurveIntersection

本项目主要用于计算和检测几何模型（如三角网格、B-spline patch 等）的相交与自交情况。

## 🚀 使用方法

在项目源码目录下，按以下三步即可完成编译与运行：

### 1. 配置项目
```bash
cmake ..
```

### 2. 编译构建
```bash
cmake --build . --config Release
```

### 3. 运行程序
编译完成后，进入 `Release` 目录，双击运行 `CurveIntersection.exe`。
根据命令行窗口的提示，输入测试用例的文件名即可执行。

---

## ⚠️ 测试文件命名规范

测试文件开头的前缀严格对应了不同的几何求交场景，输入文件名时请注意区分：

* **`tri...`** ：代表三角自交 (Triangle self-intersection)
* **`patch...`** ：代表 B-spline patch 自交 (Patch self-intersection)
* **`tritri...`** ：代表三角邻交 (Triangle-triangle adjacent intersection)
* **`patchpatch...`** ：代表 B-spline patch 邻交 (Patch-patch adjacent intersection)
* **`tripatch...`** ：代表混合邻交 (Mixed adjacent intersection)
