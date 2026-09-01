$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent (Split-Path -Parent $PSCommandPath)
python (Join-Path $Root "example/create_toy_data.py") --verify-only
python (Join-Path $Root "example/run_toy.py")
python (Join-Path $Root "scripts/check_toy_output.py") `
  --observed (Join-Path $Root "example/output/toy_run") `
  --expected (Join-Path $Root "example/expected") `
  --data (Join-Path $Root "example/data")
