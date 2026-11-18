$ErrorActionPreference = 'Stop'

node --check scripts/benchmark_gate.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_gate.mjs syntax check failed' }

node --check scripts/benchmark_score.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_score.mjs syntax check failed' }

node --check scripts/benchmark_increment.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_increment.mjs syntax check failed' }

node --check scripts/audit_gse154200.mjs
if ($LASTEXITCODE -ne 0) { throw 'audit_gse154200.mjs syntax check failed' }

node --check scripts/benchmark_gse154200.mjs
if ($LASTEXITCODE -ne 0) { throw 'benchmark_gse154200.mjs syntax check failed' }

if (Get-Command Rscript -ErrorAction SilentlyContinue) {
  Rscript -e "invisible(parse(file='scripts/generate_peerj_figures_supplements.R')); cat('R syntax check passed\n')"
  if ($LASTEXITCODE -ne 0) { throw 'generate_peerj_figures_supplements.R syntax check failed' }
} else {
  Write-Warning 'Rscript not found; R syntax check skipped.'
}

$files = Get-ChildItem -File -Recurse | Where-Object { $_.FullName -notmatch '[\\/]\.git[\\/]' -and $_.Name -ne 'SHA256_MANIFEST.csv' }
$rootPath = (Get-Location).Path
$rows = foreach ($file in $files) {
  [pscustomobject]@{
    path = $file.FullName.Substring($rootPath.Length).TrimStart('\', '/').Replace('\\', '/')
    bytes = $file.Length
    sha256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
  }
}
$rows | Sort-Object path | Export-Csv -LiteralPath 'SHA256_MANIFEST.csv' -NoTypeInformation -Encoding utf8
Write-Host "Repository verification passed. Manifest: SHA256_MANIFEST.csv"
