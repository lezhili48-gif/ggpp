# Inputs directory

请放入以下 3 个文件（文件名保持一致）：

1. `Data_Extraction_Template.xlsx`
2. `Meta_Analysis_R_Code.R`
3. `PRISMA_Protocol_Biochar_DNRA_Denitrification.pdf`

然后运行：

```bash
bash scripts/bootstrap_from_inputs.sh
```

脚本会检查文件完整性，并生成：
- `config/source_manifest.yaml`
- `config/workflow_todo.md`
