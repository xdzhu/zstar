# ZStar 的 GitHub 与 PyPI 更新流程

本文档用于以后自行发布 ZStar。先设置仓库根目录：

```powershell
$RepoRoot = "C:\path\to\zstar"
Set-Location $RepoRoot
```

## 基本约定

- `examples/` 是 GitHub 上公开的可复现实例库；它不打入 wheel，源码包是否
  包含它由构建清单决定，发布前应明确检查。
- `dist/`、`build/` 和 `*.egg-info/` 是本地构建产物，不提交 GitHub。
- 调度脚本由 ZStar CLI 按系统自动生成；仓库不再维护站点专用的旧作业模板。
- 开发和测试使用 `zstar-test` 环境。
- 构建和上传可以使用安装了 `build` 与 `twine` 的独立发布环境。
- PyPI 已发布的文件不可覆盖。

## 什么时候需要修改版本号

仅向 GitHub 提交代码或文档时，不必修改版本号。

只要希望 PyPI 页面或安装包发生变化，即使只修改 README/项目描述，也必须发布一个新版本。PyPI 不允许用同一版本号重新上传 wheel 或源码包。

建议遵循：

- 修复错误或只改文档：补丁版本，例如 `0.1.0 -> 0.1.1`。
- 增加向后兼容的新功能：次版本，例如 `0.1.0 -> 0.2.0`。
- 引入不兼容接口：主版本。

版本号需要同时修改：

- `pyproject.toml` 中的 `project.version`
- `zstar/__init__.py` 中的 `__version__`
- `CHANGELOG.md` 中的发布记录

## 1. 检查工作区

```powershell
conda activate zstar-test
Set-Location $RepoRoot

git status --short
git diff --check
```

确认没有误加入：

- `dist/`
- `build/`
- 临时输出和账号凭据

## 2. 运行本地测试

```powershell
python -m compileall -q zstar tests
python -m unittest discover -v
python -m zstar.cli --version
python -m zstar.cli --help
python -m zstar.cli workflow run --help
python -m zstar.cli agent-skill path
python -m zstar.cli agent-skill preflight --root . --lane bec --dim bulk
```

需要外部程序的例子应在仓库的 `examples/` 中验证。至少确认：

- `zstar bec post` 能处理已有极化结果。
- 默认绝缘性门控只对 `0.no-move` 执行一次普通 `--band`。
- `zstar bec stat` 能识别完成、失败和恢复状态。
- 新旧 PYATB 环境都能读取电子介电张量。

## 3. 更新中英文 README 与 PDF

编辑：

- `README.md`
- `README.zh-CN.md`
- `README_PYPI.md`

重新渲染：

```powershell
node docs\render_readme_pdfs.mjs
```

输出：

- `docs/README.en.pdf`
- `docs/README.zh-CN.pdf`

应把 PDF 转成图片进行目视检查，确认没有截断、重叠、乱码或过宽表格：

```powershell
pdftoppm -png -r 120 docs\README.en.pdf tmp\pdfs\README-en
pdftoppm -png -r 120 docs\README.zh-CN.pdf tmp\pdfs\README-zh
```

## 4. 检查 PyPI 描述

由于仓库可能是私有的，`README_PYPI.md` 不应使用：

```html
<img src="docs/logo.png">
```

也不应使用私有仓库的 `raw.githubusercontent.com` 地址。PyPI 访问者没有仓库登录权限，图片会失效。

可选方案：

1. 不在 PyPI 描述中显示 logo，这是当前默认方案。
2. 将 logo 放在独立的公开仓库或长期稳定的公共 HTTPS 静态资源服务中。
3. 把该公开 URL 写入 `README_PYPI.md`。

PyPI 不会为安装包中的 `docs/logo.png` 提供可直接嵌入项目页面的稳定资源 URL。

## 5. 清理旧构建并打包

不要删除源码目录。只清理已确认位于项目根目录下的构建产物：

```powershell
Set-Location $RepoRoot
Remove-Item -Recurse -Force -LiteralPath .\dist -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force -LiteralPath .\build -ErrorAction SilentlyContinue
Get-ChildItem -Directory -Filter *.egg-info | Remove-Item -Recurse -Force

python -m build
```

`dist/` 中应出现：

```text
zstar-X.Y.Z-py3-none-any.whl
zstar-X.Y.Z.tar.gz
```

## 6. 检查发布包内容

```powershell
python -m twine check dist\*
python -m zipfile -l dist\zstar-X.Y.Z-py3-none-any.whl
```

确认：

- 包中有 `zstar/` 模块和必要文档。
- 包中有 `zstar/agent_skills/run-zstar-workflows/`，且不含
  `__pycache__` 或 `*.pyc`。
- 包中没有 `examples/`、`dist/`、远端计算输出或凭据。
- `README_PYPI.md` 渲染检查通过。

建议创建临时环境做安装测试：

```powershell
python -m venv tmp\release-smoke
tmp\release-smoke\Scripts\python -m pip install --upgrade pip
tmp\release-smoke\Scripts\python -m pip install dist\zstar-X.Y.Z-py3-none-any.whl
tmp\release-smoke\Scripts\zstar --version
tmp\release-smoke\Scripts\zstar --help
tmp\release-smoke\Scripts\zstar skill install --dest tmp\release-smoke-skill
$codexHome = if ($env:CODEX_HOME) { $env:CODEX_HOME } else { Join-Path $HOME '.codex' }
python (Join-Path $codexHome 'skills\.system\skill-creator\scripts\quick_validate.py') tmp\release-smoke-skill\run-zstar-workflows
```

## 7. 提交并推送 GitHub

只暂存本次需要的文件：

```powershell
git add README.md README.zh-CN.md README_PYPI.md CHANGELOG.md pyproject.toml
git add zstar tests docs README*.md MANIFEST.in
git status --short
git diff --cached --check
```

提交并推送：

```powershell
git commit -m "Release ZStar X.Y.Z"
git push origin main
```

`examples/` 是公开案例内容，应在提交前检查其 manifest、输入路径和参考结果；
`dist/` 仍不提交。

## 8. 上传 PyPI

推荐使用 PyPI API token，不要把 token 写入仓库。

可先测试：

```powershell
python -m twine upload --repository testpypi dist\*
```

正式上传：

```powershell
python -m twine upload dist\*
```

用户名使用：

```text
__token__
```

密码使用 PyPI 生成的 token。

## 9. 发布后验证

等待 PyPI 页面刷新后：

```powershell
python -m venv tmp\pypi-smoke
tmp\pypi-smoke\Scripts\python -m pip install --upgrade pip
tmp\pypi-smoke\Scripts\python -m pip install --no-cache-dir zstar==X.Y.Z
tmp\pypi-smoke\Scripts\zstar --version
tmp\pypi-smoke\Scripts\zstar --help
```

同时检查：

- PyPI 项目描述中的表格、代码块和链接。
- GitHub 中英 README 与 logo。
- GitHub PDF 链接。
- 仓库提交中包含整理后的 `examples/`，但不包含 `dist/`、scratch 输出或凭据。

## 常见问题

### PyPI 提示文件已经存在

同一版本号不能重复上传。提高版本号、重新构建并再次上传。

### GitHub 能显示 logo，但 PyPI 不能

相对路径图片只对仓库页面有效。私有 GitHub 原始文件对未登录的 PyPI 访问者不可见，应使用公开 HTTPS URL 或不显示 logo。

### 只改 README，是否需要更新软件版本

只推送 GitHub：不需要。

要让 PyPI 页面同步新的 README：需要，因为必须上传一个新的、版本号不同的发布包。
