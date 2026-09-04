import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath, pathToFileURL } from "node:url";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const runtimeModules = path.join(
  os.homedir(),
  ".cache",
  "codex-runtimes",
  "codex-primary-runtime",
  "dependencies",
  "node",
  "node_modules",
);
const modulePaths = [
  path.join(repoRoot, "node_modules"),
  process.env.CODEX_NODE_MODULES,
  runtimeModules,
].filter(Boolean);

async function loadPackage(name) {
  const resolved = require.resolve(name, { paths: modulePaths });
  return import(pathToFileURL(resolved).href);
}

const { marked } = await loadPackage("marked");
const css = fs.readFileSync(path.join(repoRoot, "docs", "readme-print.css"), "utf8");
const pyproject = fs.readFileSync(path.join(repoRoot, "pyproject.toml"), "utf8");
const versionMatch = pyproject.match(/^version\s*=\s*"([^"]+)"/m);
if (!versionMatch) {
  throw new Error("Could not read the project version from pyproject.toml.");
}
const tempDir = path.join(repoRoot, "tmp", "pdfs", "readme-render");
fs.mkdirSync(tempDir, { recursive: true });

const documents = [
  ["README.md", "docs/README.en.pdf", "ZStar English Manual"],
  ["README.zh-CN.md", "docs/README.zh-CN.pdf", "ZStar 中文手册"],
];

const browserCandidates = [
  process.env.CHROME_PATH,
  "C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe",
  "C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe",
  "C:\\Program Files\\Microsoft\\Edge\\Application\\msedge.exe",
  "google-chrome",
  "chromium",
].filter(Boolean);
const browserExecutable = browserCandidates.find((candidate) => {
  if (path.isAbsolute(candidate)) {
    return fs.existsSync(candidate);
  }
  return spawnSync(candidate, ["--version"], { stdio: "ignore" }).status === 0;
});
if (!browserExecutable) {
  throw new Error("Chrome, Edge, or Chromium was not found. Set CHROME_PATH.");
}

for (const [sourceName, outputName, title] of documents) {
  const markdown = fs
    .readFileSync(path.join(repoRoot, sourceName), "utf8")
    .replace(
      /<img alt="PyPI"[^>]*>/g,
      `<span class="pdf-badge badge-pypi">ZStar v${versionMatch[1]}</span>`,
    )
    .replace(
      /<img alt="Python"[^>]*>/g,
      '<span class="pdf-badge badge-python">Python 3.9+</span>',
    )
    .replace(
      /<img alt="(?:License|许可证)"[^>]*>/g,
      '<span class="pdf-badge badge-license">License GPL-3.0</span>',
    );
  const content = marked.parse(markdown, { gfm: true });
  const html = `<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <base href="${pathToFileURL(`${repoRoot}${path.sep}`).href}">
  <title>${title}</title>
  <style>${css}</style>
</head>
<body><main>${content}</main></body>
</html>`;
  const htmlPath = path.join(tempDir, `${path.basename(sourceName)}.html`);
  fs.writeFileSync(htmlPath, html, "utf8");

  const outputPath = path.join(repoRoot, outputName);
  const profilePath = path.join(tempDir, `chrome-${path.basename(sourceName)}`);
  const result = spawnSync(browserExecutable, [
    "--headless=new",
    "--disable-gpu",
    "--disable-extensions",
    "--disable-background-networking",
    "--no-pdf-header-footer",
    "--virtual-time-budget=1500",
    `--user-data-dir=${profilePath}`,
    `--print-to-pdf=${outputPath}`,
    pathToFileURL(htmlPath).href,
  ], {
    encoding: "utf8",
    timeout: 120000,
  });
  if (result.status !== 0 || !fs.existsSync(outputPath)) {
    throw new Error(
      `PDF rendering failed for ${sourceName}: ${result.stderr || result.stdout}`
    );
  }
  process.stdout.write(`Rendered ${outputName}\n`);
}
