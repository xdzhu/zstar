# How to update `zstar` on PyPI

约定：

- 开发和本地验证环境：`zstar-test` 或当前开发环境
- 打包和上传环境：`base`
- 版本号格式：`X.Y.Z`，例如 `0.0.8`
- `examples/` 只用于本地验证，不提交到 GitHub
- `dist/` 和 `build/` 是构建产物，不提交到 GitHub

## 1. 更新代码和版本号

```powershell
conda activate zstar-test
cd D:\Work\Code\zstar
```

确认本地例子和核心命令可以正常运行后，修改两处版本号：

- `pyproject.toml`

```toml
[project]
name = "zstar"
version = "X.Y.Z"
```

- `zstar/__init__.py`

```python
__version__ = "X.Y.Z"
```

同时建议在 `CHANGELOG.md` 增加这一版的变更说明。

## 2. 本地验证

```powershell
python -m zstar.cli --version
python -m zstar.cli --help
python -m compileall zstar
```

如需验证示例，请在 `examples/` 目录中运行对应命令。该目录已加入 `.gitignore`，只保留在本地。

## 3. 构建发布包

回到 `base` 环境：

```powershell
conda deactivate
cd D:\Work\Code\zstar

Remove-Item -Recurse -Force dist, build -ErrorAction SilentlyContinue
python -m build
```

构建成功后，`dist/` 中应包含：

- `zstar-X.Y.Z.tar.gz`
- `zstar-X.Y.Z-py3-none-any.whl`

## 4. 检查发布包

```powershell
python -m twine check dist/*
```

可选：先上传到 TestPyPI。

```powershell
python -m twine upload --repository testpypi dist/*.whl
```

## 5. 上传到正式 PyPI

```powershell
python -m twine upload dist/*
```

上传成功后验证：

```powershell
pip install -U zstar
zstar --version
```

确认输出版本为 `X.Y.Z`。

## 6. 同步到 GitHub

发布前确认不要提交本地产物：

```powershell
git status --short
```

`examples/`、`dist/`、`build/` 不应出现在待提交文件中。

提交并推送：

```powershell
git add .
git commit -m "Release zstar X.Y.Z"
git push origin main
```
