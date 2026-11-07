# RMSDAnalyzer

동일 단백질의 두 PDB 구조를 비교하고, 잔기별 Cα RMSD와 backbone RMSD를
계산하는 도구입니다. 두 개의 PDB 파일 또는 하나의 PDB 파일에 포함된 두
`MODEL`을 입력으로 사용할 수 있습니다.

## 실행

```bash
python3 -m pip install -r requirements.txt
./run_analysis.sh
```

실행 전 [`run_analysis.sh`](run_analysis.sh)의 `USER SETTINGS` 영역에서 입력과
옵션을 수정합니다.

## Input

- 두 PDB 파일 또는 하나의 PDB 파일에 포함된 두 `MODEL`
- `align`: 서열 정렬 후 구조 정렬
- `prealigned`: 현재 좌표를 그대로 비교
- 지정 잔기 선택 가능 (`A:42`, `A:42A`)

## Output

- 전체 Cα RMSD
- chain 대응 결과
- 잔기별 Cα RMSD와 backbone RMSD
- backbone RMSD TOP 10
- 지정 잔기 결과
- `output/residue_rmsd.csv`
