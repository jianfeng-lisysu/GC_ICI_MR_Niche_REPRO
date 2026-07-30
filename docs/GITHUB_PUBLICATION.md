# Publishing this repository to GitHub

## Required public URL

`https://github.com/jianfeng-lisysu/GC_ICI_MR_Niche_REPRO`

The manuscript states that the analysis code is publicly available at this URL. Do not submit the manuscript until this URL opens for a signed-out user and the repository is confirmed as **Public**.

## Recommended Windows workflow

1. Extract the ZIP file.
2. Open the extracted `GC_ICI_MR_Niche_REPRO` folder.
3. Double-click `PUBLISH_TO_GITHUB.cmd`.
4. Sign in through the browser when GitHub CLI requests authentication.
5. The script will:
   - verify the manifest and SHA-256 checksums;
   - initialize the local Git repository;
   - create the repository if it does not exist;
   - change it to Public if it already exists as Private;
   - push the `main` branch;
   - verify public anonymous access.
6. After completion, double-click `VERIFY_PUBLIC_REPOSITORY.cmd`.
7. Open the repository in a private/incognito browser window and confirm that `README.md` and `scripts/` are visible.

## Required command-line tools

When Git or GitHub CLI is missing, install them in PowerShell:

```powershell
winget install --id Git.Git -e
winget install --id GitHub.cli -e
```

Close and reopen PowerShell after installation.

## Manual equivalent

```powershell
gh auth login --web --git-protocol https
gh auth setup-git

git init
git branch -M main
git add -A
git commit -m "Initial public reproducibility release"

gh repo create jianfeng-lisysu/GC_ICI_MR_Niche_REPRO `
  --public `
  --source . `
  --remote origin `
  --push
```

For an existing private repository:

```powershell
gh repo edit jianfeng-lisysu/GC_ICI_MR_Niche_REPRO `
  --visibility public `
  --accept-visibility-change-consequences

git remote set-url origin `
  https://github.com/jianfeng-lisysu/GC_ICI_MR_Niche_REPRO.git

git push -u origin main
```

## Submission gate

The repository is ready for citation only after all of the following are true:

- the anonymous GitHub API returns the repository;
- `private` is `false`;
- the public `README.md` opens;
- the `scripts` directory is visible;
- the default branch contains the latest cleaned repository files.
