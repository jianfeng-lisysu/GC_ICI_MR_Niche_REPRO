$ErrorActionPreference = "Stop"

$Owner = "jianfeng-lisysu"
$Repository = "GC_ICI_MR_Niche_REPRO"
$FullName = "$Owner/$Repository"
$RepositoryUrl = "https://github.com/$FullName"

Write-Host "============================================================"
Write-Host "Anonymous public-repository verification"
Write-Host "============================================================"
Write-Host "Repository: $RepositoryUrl"
Write-Host ""

try {
    $Metadata = Invoke-RestMethod `
        -Uri "https://api.github.com/repos/$FullName" `
        -Headers @{
            "Accept" = "application/vnd.github+json"
            "User-Agent" = "GC-ICI-MR-Niche-publication-check"
        }

    if ($Metadata.private -ne $false) {
        throw "GitHub API did not report private=false."
    }

    if ($Metadata.full_name -ne $FullName) {
        throw "GitHub API returned '$($Metadata.full_name)' instead of '$FullName'."
    }

    $Readme = Invoke-WebRequest `
        -Uri "https://raw.githubusercontent.com/$FullName/main/README.md" `
        -Headers @{
            "User-Agent" = "GC-ICI-MR-Niche-publication-check"
        } `
        -UseBasicParsing

    if ($Readme.StatusCode -ne 200) {
        throw "The public README returned HTTP $($Readme.StatusCode)."
    }

    if ($Readme.Content -notmatch "GC_ICI_MR_Niche_REPRO") {
        throw "The public README does not contain the expected repository name."
    }

    Write-Host "PASS" -ForegroundColor Green
    Write-Host "Visibility: Public"
    Write-Host "Anonymous API access: OK"
    Write-Host "Public README: OK"
    Write-Host "URL: $RepositoryUrl"
    Start-Process $RepositoryUrl
}
catch {
    Write-Host "FAIL" -ForegroundColor Red
    Write-Host $_.Exception.Message
    Write-Host ""
    Write-Host "Do not submit the manuscript while this check fails."
    exit 1
}

Read-Host "Press Enter to close"
