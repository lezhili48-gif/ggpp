#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="${1:-inputs}"

required_files=(
  "Data_Extraction_Template.xlsx"
  "Meta_Analysis_R_Code.R"
  "PRISMA_Protocol_Biochar_DNRA_Denitrification.pdf"
)

missing=0
for f in "${required_files[@]}"; do
  if [[ ! -f "${INPUT_DIR}/${f}" ]]; then
    echo "[MISSING] ${INPUT_DIR}/${f}"
    missing=1
  else
    echo "[OK] ${INPUT_DIR}/${f}"
  fi
done

if [[ ${missing} -ne 0 ]]; then
  echo "\n请先把以上缺失文件放到 ${INPUT_DIR}/ 后再运行。" >&2
  exit 1
fi

mkdir -p artifacts logs config

cat > config/source_manifest.yaml <<MANIFEST
repository: ggpp
branch: main
inputs:
  extraction_template: "${INPUT_DIR}/Data_Extraction_Template.xlsx"
  analysis_code: "${INPUT_DIR}/Meta_Analysis_R_Code.R"
  prisma_protocol: "${INPUT_DIR}/PRISMA_Protocol_Biochar_DNRA_Denitrification.pdf"
MANIFEST

cat > config/workflow_todo.md <<'TODO'
# Workflow TODO (from provided inputs)

## 1) Protocol parsing
- Read PRISMA protocol and extract:
  - PICO/PECO definition
  - Inclusion/exclusion criteria
  - Primary/secondary outcomes
  - Planned subgroup/sensitivity analysis

## 2) Data schema mapping
- From Data_Extraction_Template.xlsx, map each sheet/column to:
  - field name
  - type/unit
  - missing-value policy
  - evidence trace field (page/table/quote)

## 3) Statistical reproducibility
- Validate Meta_Analysis_R_Code.R:
  - package versions
  - fixed random seed
  - I/O paths standardized
  - outputs (forest/funnel/tables) to artifacts/

## 4) Manuscript binding
- Link analysis outputs to manuscript placeholders
- Generate result snippets with figure/table references
TODO

echo "\n生成完成："
echo "- config/source_manifest.yaml"
echo "- config/workflow_todo.md"
echo "\n下一步："
echo "1) bash scripts/bootstrap_from_inputs.sh"
echo "2) 把 R 脚本输出目录改为 artifacts/"
