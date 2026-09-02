# Workspace Sync

Primary local workspace:

`E:\AFDM\afdm_thesis_workspace`

Primary branch:

`codex/thesis-workspace`

Remote:

`git@github.com:jiaqishi-cmd/AFDM_pre_chirp.git`

## Daily Workflow

Before leaving one computer:

```bash
git status
git add <files-to-sync>
git commit -m "Work checkpoint"
git push
```

After opening the other computer:

```bash
git checkout codex/thesis-workspace
git pull
```

## What Belongs In Git

- MATLAB source code and experiment scripts
- Research skills
- Notes, idea logs, reading summaries, and experiment plans
- Small configuration or metadata files

## What Usually Stays Out

- Large `.mat` files
- Generated figures
- Logs
- Temporary simulation results
- Local editor settings

