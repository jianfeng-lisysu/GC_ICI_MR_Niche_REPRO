$ErrorActionPreference = "Stop"

$Owner = "jianfeng-lisysu"
$Repository = "GC_ICI_MR_Niche_REPRO"
$FullName = "$Owner/$Repository"
$RepositoryUrl = "https://github.com/$FullName"
$RepositoryDescription = "Reproducibility code for the FCN1/TNFSF12 stromal-myeloid gastric cancer study."

Set-Location $PSScriptRoot

function Stop-WithMessage {
    param([string]$Message)
    Write-Host ""
    Write-Host $Message -ForegroundColor Red
    Write-Host ""
    Read-Host "Press Enter to close"
    exit 1
}

Write-Host "============================================================"
Write-Host "Publish GC_ICI_MR_Niche_REPRO to GitHub"
Write-Host "============================================================"
Write-Host "Target: $RepositoryUrl"
Write-Host ""

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    Stop-WithMessage @"
Git is not installed or is not available in PATH.

Open PowerShell as a normal user and run:
winget install --id Git.Git -e

Then close and reopen PowerShell before running this script again.
"@
}

if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    Stop-WithMessage @"
GitHub CLI is not installed or is not available in PATH.

Open PowerShell as a normal user and run:
winget install --id GitHub.cli -e

Then close and reopen PowerShell before running this script again.
"@
}

Write-Host "Checking GitHub authentication..."
gh auth status 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "A browser window will open for GitHub authentication."
    gh auth login --web --git-protocol https
    if ($LASTEXITCODE -ne 0) {
        Stop-WithMessage "GitHub authentication failed."
    }
}

gh auth setup-git
if ($LASTEXITCODE -ne 0) {
    Stop-WithMessage "Could not configure Git to use GitHub CLI credentials."
}

$AuthenticatedUser = gh api user --jq ".login"
if ($LASTEXITCODE -ne 0) {
    Stop-WithMessage "Could not identify the authenticated GitHub account."
}

if ($AuthenticatedUser.Trim().ToLowerInvariant() -ne $Owner.ToLowerInvariant()) {
    Stop-WithMessage @"
The active GitHub account is '$AuthenticatedUser', but this repository must be published under '$Owner'.

Run:
gh auth switch --user $Owner

Then run this script again.
"@
}

# Verify that the package is internally complete before publication.
$ManifestPath = Join-Path $PSScriptRoot "FILE_MANIFEST.txt"
$ChecksumPath = Join-Path $PSScriptRoot "CHECKSUMS_SHA256.txt"

if (-not (Test-Path $ManifestPath)) {
    Stop-WithMessage "FILE_MANIFEST.txt is missing."
}
if (-not (Test-Path $ChecksumPath)) {
    Stop-WithMessage "CHECKSUMS_SHA256.txt is missing."
}

$MissingManifestFiles = @()
Get-Content $ManifestPath | ForEach-Object {
    $RelativePath = $_.Trim()
    if ($RelativePath -and -not (Test-Path (Join-Path $PSScriptRoot $RelativePath))) {
        $MissingManifestFiles += $RelativePath
    }
}
if ($MissingManifestFiles.Count -gt 0) {
    Stop-WithMessage ("Files listed in FILE_MANIFEST.txt are missing:`n" + ($MissingManifestFiles -join "`n"))
}

$ChecksumFailures = @()
Get-Content $ChecksumPath | ForEach-Object {
    if ($_ -match "^([0-9a-fA-F]{64})  (.+)$") {
        $Expected = $Matches[1].ToLowerInvariant()
        $RelativePath = $Matches[2]
        $FullPath = Join-Path $PSScriptRoot $RelativePath

        if (-not (Test-Path $FullPath)) {
            $ChecksumFailures += "$RelativePath (missing)"
        }
        else {
            $Observed = (Get-FileHash -Path $FullPath -Algorithm SHA256).Hash.ToLowerInvariant()
            if ($Observed -ne $Expected) {
                $ChecksumFailures += "$RelativePath (checksum mismatch)"
            }
        }
    }
}
if ($ChecksumFailures.Count -gt 0) {
    Stop-WithMessage ("Checksum validation failed:`n" + ($ChecksumFailures -join "`n"))
}

if (-not (Test-Path (Join-Path $PSScriptRoot ".git"))) {
    git init
    if ($LASTEXITCODE -ne 0) {
        Stop-WithMessage "git init failed."
    }
}

git branch -M main
if ($LASTEXITCODE -ne 0) {
    Stop-WithMessage "Could not set the local branch to main."
}

$GitName = git config user.name
if (-not $GitName) {
    $GitName = Read-Host "Enter the Git commit author name"
    git config user.name "$GitName"
}

$GitEmail = git config user.email
if (-not $GitEmail) {
    $GitEmail = Read-Host "Enter the Git commit author email"
    git config user.email "$GitEmail"
}

git add -A
if ($LASTEXITCODE -ne 0) {
    Stop-WithMessage "git add failed."
}

$PendingChanges = git status --porcelain
if ($PendingChanges) {
    git commit -m "Initial public reproducibility release"
    if ($LASTEXITCODE -ne 0) {
        Stop-WithMessage "git commit failed."
    }
}
else {
    Write-Host "No uncommitted changes were found."
}

$RepositoryExists = $false
gh repo view $FullName --json nameWithOwner,visibility,url 1>$null 2>$null
if ($LASTEXITCODE -eq 0) {
    $RepositoryExists = $true
}

if (-not $RepositoryExists) {
    Write-Host "Creating the public GitHub repository..."
    gh repo create $FullName `
        --public `
        --description $RepositoryDescription `
        --source . `
        --remote origin `
        --push

    if ($LASTEXITCODE -ne 0) {
        Stop-WithMessage "GitHub repository creation or initial push failed."
    }
}
else {
    Write-Host "The GitHub repository already exists."

    $Visibility = gh repo view $FullName --json visibility --jq ".visibility"
    if ($LASTEXITCODE -ne 0) {
        Stop-WithMessage "Could not read repository visibility."
    }

    if ($Visibility.Trim().ToUpperInvariant() -ne "PUBLIC") {
        Write-Host "Changing repository visibility to Public..."
        gh repo edit $FullName `
            --visibility public `
            --accept-visibility-change-consequences

        if ($LASTEXITCODE -ne 0) {
            Stop-WithMessage "Could not change repository visibility to Public."
        }
    }

    $OriginUrl = "https://github.com/$FullName.git"
    $ExistingOrigin = git remote get-url origin 2>$null

    if ($LASTEXITCODE -ne 0) {
        git remote add origin $OriginUrl
    }
    elseif ($ExistingOrigin.Trim() -ne $OriginUrl) {
        git remote set-url origin $OriginUrl
    }

    git push -u origin main
    if ($LASTEXITCODE -ne 0) {
        Stop-WithMessage @"
Push failed. The remote repository may contain commits that are not present locally.
Do not use a force push without reviewing the remote history.
"@
    }
}

# Set stable metadata after the push.
gh repo edit $FullName `
    --description $RepositoryDescription `
    --add-topic gastric-cancer `
    --add-topic mendelian-randomization `
    --add-topic spatial-transcriptomics `
    --add-topic immunotherapy `
    --add-topic reproducibility

if ($LASTEXITCODE -ne 0) {
    Write-Host "Warning: repository topics or description could not be updated." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "Running authenticated visibility check..."
$Visibility = gh repo view $FullName --json visibility --jq ".visibility"
if ($Visibility.Trim().ToUpperInvariant() -ne "PUBLIC") {
    Stop-WithMessage "The repository exists, but its visibility is not Public."
}

Write-Host "Running anonymous GitHub API check..."
try {
    $PublicMetadata = Invoke-RestMethod `
        -Uri "https://api.github.com/repos/$FullName" `
        -Headers @{
            "Accept" = "application/vnd.github+json"
            "User-Agent" = "GC-ICI-MR-Niche-publication-check"
        }

    if ($PublicMetadata.private -ne $false) {
        Stop-WithMessage "Anonymous API check did not confirm a public repository."
    }

    if ($PublicMetadata.full_name -ne $FullName) {
        Stop-WithMessage "Anonymous API check returned the wrong repository."
    }
}
catch {
    Stop-WithMessage ("Anonymous GitHub API check failed: " + $_.Exception.Message)
}

Write-Host ""
Write-Host "SUCCESS" -ForegroundColor Green
Write-Host "Repository: $RepositoryUrl"
Write-Host "Visibility: Public"
Write-Host "Default branch: main"
Write-Host ""
Write-Host "Open this URL in a signed-out/private browser window before submission."
Start-Process $RepositoryUrl
Read-Host "Press Enter to close"
