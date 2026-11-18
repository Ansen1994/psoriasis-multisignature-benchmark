$ErrorActionPreference = 'Stop'
$repoRoot = $PSScriptRoot
Push-Location $repoRoot
try {

if (-not (Get-Command node -ErrorAction SilentlyContinue)) {
  throw 'Node.js is required to run the analysis workflow.'
}

$requiredInputs = @(
  'sources/raw/GPL570/GPL570.annot.gz',
  'sources/raw/GPL96/GPL96_full.txt',
  'sources/raw/GSE13355/GSE13355_series_matrix.txt.gz',
  'sources/raw/GSE14905/GSE14905_series_matrix.txt.gz',
  'sources/raw/GSE154200/GSE154200_series_matrix.txt.gz',
  'sources/raw/GSE154200/GSE154200_countMat.txt.gz',
  'sources/raw/GSE154200/GSE154200_samples_brief.txt',
  'sources/raw/GSE212476/GSE212476_submission_fpkm.csv.gz',
  'sources/raw/GSE30999/GSE30999_series_matrix.txt.gz',
  'sources/raw/GSE41664/GSE41664_series_matrix.txt.gz',
  'sources/raw/GSE6710/GSE6710_series_matrix.txt.gz',
  'sources/raw/GSE78097/GSE78097_series_matrix.txt.gz',
  'sources/raw/msigdb/c5.go.v2025.1.Hs.symbols.gmt',
  'sources/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt'
)

$missing = $requiredInputs | Where-Object { -not (Test-Path -LiteralPath $_) }
if ($missing) {
  throw "Missing public inputs:`n$($missing -join "`n")"
}

node scripts/benchmark_gate.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_gate.mjs failed' }

node scripts/benchmark_score.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_score.mjs failed' }

node scripts/benchmark_increment.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_increment.mjs failed' }

node scripts/audit_gse154200.mjs
if ($LASTEXITCODE -ne 0) { throw 'audit_gse154200.mjs failed' }

node scripts/benchmark_gse154200.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_gse154200.mjs failed' }

if (Get-Command Rscript -ErrorAction SilentlyContinue) {
  Rscript scripts/generate_peerj_figures_supplements.R
  if ($LASTEXITCODE -ne 0) { throw 'generate_peerj_figures_supplements.R failed' }
} else {
  Write-Warning 'Rscript not found; analysis outputs were generated, but figures and supplements were not rebuilt.'
}

} finally {
  Pop-Location
}
