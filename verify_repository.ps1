param(
    [switch]$UpdateManifest
)

$ErrorActionPreference = 'Stop'
$repoRoot = $PSScriptRoot
$manifestPath = Join-Path $repoRoot 'SHA256_MANIFEST.csv'
Push-Location $repoRoot
try {
    if (-not (Get-Command node -ErrorAction SilentlyContinue)) {
        throw 'Node.js is required for the JavaScript syntax checks.'
    }

    $nodeScripts = @(
        'scripts/benchmark_gate.mjs',
        'scripts/benchmark_score.mjs',
        'scripts/benchmark_increment.mjs',
        'scripts/audit_gse154200.mjs',
        'scripts/benchmark_gse154200.mjs'
    )
    foreach ($script in $nodeScripts) {
        node --check $script
        if ($LASTEXITCODE -ne 0) { throw "$script syntax check failed" }
    }

    if (Get-Command Rscript -ErrorAction SilentlyContinue) {
        Rscript -e "invisible(parse(file='scripts/generate_peerj_figures_supplements.R')); cat('R syntax check passed\n')"
        if ($LASTEXITCODE -ne 0) { throw 'generate_peerj_figures_supplements.R syntax check failed' }
    } else {
        Write-Warning 'Rscript not found; R syntax check skipped.'
    }

    if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
        throw 'Git is required to determine the repository file set.'
    }
    $gitPaths = git ls-files --cached --others --exclude-standard
    if ($LASTEXITCODE -ne 0) { throw 'git ls-files failed.' }

    $rows = foreach ($relativePath in $gitPaths) {
        $normalized = $relativePath.Replace('\', '/')
        if ($normalized -eq 'SHA256_MANIFEST.csv') { continue }
        $fullPath = Join-Path $repoRoot $relativePath
        if (-not (Test-Path -LiteralPath $fullPath -PathType Leaf)) {
            throw "Tracked file is missing from disk: $normalized"
        }
        [pscustomobject]@{
            path = $normalized
            bytes = (Get-Item -LiteralPath $fullPath).Length
            sha256 = (Get-FileHash -LiteralPath $fullPath -Algorithm SHA256).Hash.ToLowerInvariant()
        }
    }
    $rows = @($rows | Sort-Object path)

    if ($UpdateManifest) {
        $csv = (($rows | ConvertTo-Csv -NoTypeInformation) -join "`r`n") + "`r`n"
        [System.IO.File]::WriteAllText($manifestPath, $csv, (New-Object System.Text.UTF8Encoding($false)))
        Write-Host "Manifest updated: $manifestPath"
        return
    }

    if (-not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) {
        throw 'SHA256_MANIFEST.csv is missing; run with -UpdateManifest after reviewing changes.'
    }

    $expected = @(Import-Csv -LiteralPath $manifestPath)
    $currentByPath = @{}
    foreach ($row in $rows) { $currentByPath[$row.path] = $row }
    $expectedByPath = @{}
    foreach ($row in $expected) { $expectedByPath[$row.path] = $row }

    $allPaths = @($currentByPath.Keys + $expectedByPath.Keys | Sort-Object -Unique)
    $differences = foreach ($path in $allPaths) {
        if (-not $currentByPath.ContainsKey($path)) {
            "missing from repository: $path"
            continue
        }
        if (-not $expectedByPath.ContainsKey($path)) {
            "not listed in manifest: $path"
            continue
        }
        $actual = $currentByPath[$path]
        $recorded = $expectedByPath[$path]
        if ($actual.bytes -ne $recorded.bytes -or $actual.sha256 -ne $recorded.sha256) {
            "changed: $path"
        }
    }

    if ($differences) {
        throw "Repository manifest mismatch:`n$($differences -join "`n")`nRun with -UpdateManifest only after reviewing the changes."
    }

    Write-Host "Repository verification passed. $($rows.Count) files match SHA256_MANIFEST.csv."
} finally {
    Pop-Location
}
