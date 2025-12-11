# How to update `zstar` on PyPI（超简版）

> 约定：
>
> * **开发环境**：`zstar-dev`（平时调代码用）
> * **打包 / 上传环境**：`base`

版本号用占位 `X.Y.Z`，比如 `0.0.3`、`0.1.0`。

---

## 1. 在 `zstar-dev` 环境里改代码 & 版本号

```bash
conda activate zstar-test
cd D:\Work\Code\zstar
```

1. 改代码，确保本地算例里 `zstar` 能正常跑。

2. 修改版本号（两处）：

   * `pyproject.toml`：

     ```toml
     [project]
     name = "zstar"
     version = "X.Y.Z"
     ```

   * `zstar/__init__.py`：

     ```python
     __version__ = "X.Y.Z"
     ```

---

## 2. 在 `base` 环境打包

```bash
conda deactivate          # 退回 base
cd D:\Work\Code\zstar

# 清理旧构建
rm -r -fo dist,build 2>$null

# 构建 sdist + wheel
python -m build
```

构建成功后，`dist/` 里应有：

* `zstar-X.Y.Z.tar.gz`
* `zstar-X.Y.Z-py3-none-any.whl`

---


## 5. 上传到正式 PyPI

回到 `base`：

```bash
conda deactivate
cd D:\Work\Code\zstar

# python -m twine upload dist/*
# 如果只想上传 wheel，就：
python -m twine upload dist/*.whl
```

上传成功后，全局即可：

```bash
pip install -U zstar
zstar --version   # 应该显示 X.Y.Z
```

---


## 3. 上传到 TestPyPI（可选）

```bash
python -m twine upload --repository testpypi dist/*.whl
```

---

## 4. 在测试环境里验证一下（可选）

```bash
conda activate zstar-test

pip uninstall -y zstar
pip install zstar

zstar --version
```

确认输出版本是 `X.Y.Z`，核心子命令能正常 `--help`。

---

就这些。
以后发版按这个顺序：**zstar-dev 改版本 → base 打包 → testpypi（可选）→ pypi 上传**。
