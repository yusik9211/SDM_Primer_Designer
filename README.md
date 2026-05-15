# SDM Primer Designer

Site-Directed Mutagenesis(SDM)용 Double Primer 자동 설계 도구.  
Template DNA 서열과 원하는 변이를 입력하면 NEB 방식의 SDM 설계 규칙을 모두 만족하는 Forward/Reverse primer 쌍을 자동으로 디자인합니다.

---

## 파일 구성

| 파일 | 설명 |
|------|------|
| `sdm_primer_designer_v3.html` | **웹 버전 (권장)** — 브라우저에서 바로 실행, Excel/TXT 다운로드 |
| `sdm_primer_designer_v3.py`  | **Python 버전** — 터미널에서 실행, Python 3 필요 |
| `sdm_primer_designer_v2.html` | 이전 버전 (gene name 기능 없음) |
| `sdm_primer_designer_v2.py`  | 이전 버전 (gene name 기능 없음) |

---

## 빠른 시작

### HTML 버전
별도 설치 없이 `sdm_primer_designer_v3.html` 파일을 **브라우저에서 열면** 바로 사용 가능합니다.  
(인터넷 연결 필요 — SheetJS CDN을 통해 Excel 내보내기)

### Python 버전
```bash
python sdm_primer_designer_v3.py
```
Python 3.6 이상, 외부 패키지 불필요.

---

## 사용 방법 (HTML 버전)

### Step 1 — Template DNA 서열 입력
- Template DNA 서열을 5'→3' 방향으로 붙여넣기 (ATGCN만 허용)
- **Gene name** 입력 (예: `MaEgtB`) — 입력 시 primer 이름에 자동 반영
- CDS start position 설정 (AA 변이 모드 사용 시, default = 1)

### Step 2 — 변이 입력
변이 타입을 선택해 추가합니다. 여러 개를 동시에 입력할 수 있습니다.

#### DNA 변이 (DNA Mutation)
| 필드 | 설명 | 예시 |
|------|------|------|
| Position | 변이 시작 위치 (1-based) | `298` |
| Original | 원래 서열 (`-` = 삽입) | `GGC` |
| New | 새 서열 (`-` = 결실) | `GCC` |

#### 아미노산 변이 (AA Mutation)
| 필드 | 설명 | 예시 |
|------|------|------|
| AA Position | 아미노산 위치 (1-based) | `100` |
| Orig AA | 원래 아미노산 (1-letter) | `G` |
| New AA | 새 아미노산 (1-letter) | `A` |

> AA 변이 모드는 E. coli K-12 MG1655 codon usage 빈도 기준으로 최적 코돈을 자동 선택합니다.

### Step 3 — Design Primers 클릭
- **Individual Primer Pairs**: 변이마다 독립적인 primer 쌍 설계
- **Combined Primer Pair**: 변이가 2개 이상일 때 모든 변이를 하나의 primer 쌍에 포함

### 결과 확인 및 다운로드
- 각 primer 위에 이름 태그(예: `MaEgtB_G100A_F`) 표시
- Validation Checklist로 설계 규칙 충족 여부 확인
- **Download TXT**: 전체 결과를 텍스트 파일로 저장
- **Download Excel**: primer 이름, 서열, Tm, %GC 등을 Excel 파일로 저장

---

## 사용 방법 (Python 버전)

```
python sdm_primer_designer_v3.py
```

실행 후 프롬프트에 따라 순서대로 입력합니다.

```
Template DNA sequence: ATGCGT...    ← DNA 서열 붙여넣기
Gene name []:          MaEgtB       ← 유전자 이름 (Enter로 건너뛰기 가능)
CDS start [1]:         1            ← CDS 시작 위치

Design mode:
  [1] Single mutation
  [2] Multiple mutations
Choose [1/2]: 2
```

#### 변이 입력 형식
```
# DNA 변이
<position> <original> <new>
73 TAC GCG          → nt73 TAC→GCG (치환)
45 - GAATTC         → nt45에 GAATTC 삽입
30 TGCA -           → nt30부터 TGCA 결실

# DNA 변이 (nt 접두사 사용 가능)
nt 73 TAC GCG

# 아미노산 변이
aa <aa_position> <orig_AA> <new_AA>
aa 100 G A          → 100번째 Gly → Ala
aa 25 Y A           → 25번째 Tyr → Ala
```

---

## Primer 명명 규칙

```
{Gene}_{변이}_{방향}
```

| 변이 타입 | 예시 입력 | Forward 이름 | Reverse 이름 |
|----------|----------|-------------|-------------|
| AA 치환 | aa 100 G A | `MaEgtB_G100A_F` | `MaEgtB_G100A_R` |
| 단일 nt 치환 | 300 T A | `MaEgtB_T300A_F` | `MaEgtB_T300A_R` |
| 다중 nt 치환 | 80 ATG CTG | `MaEgtB_nt80ATGtoCTG_F` | `MaEgtB_nt80ATGtoCTG_R` |
| 삽입 | 45 - GAATTC | `MaEgtB_nt45insGAATTC_F` | `MaEgtB_nt45insGAATTC_R` |
| 결실 | 30 TGCA - | `MaEgtB_nt30delTGCA_F` | `MaEgtB_nt30delTGCA_R` |
| 다중 변이 (combined) | G100A + L200P | `MaEgtB_G100A_L200P_F` | `MaEgtB_G100A_L200P_R` |

Gene name 미입력 시 `G100A_F` 형식으로 표기됩니다.

---

## Primer 설계 규칙

모든 설계 결과는 아래 7가지 조건을 동시에 만족합니다.

| 규칙 | 조건 |
|------|------|
| 변이 포함 | Forward/Reverse primer 모두에 변이 포함 |
| 5' 여유 | 변이 위치가 5' 말단으로부터 **≥ 4 nt** |
| 3' 여유 | 변이 위치가 3' 말단으로부터 **≥ 8 nt** |
| Non-overlapping 3' 연장 | 각 primer의 3' 말단 중 상대 primer와 겹치지 않는 부분 **≥ 8 nt** |
| 3' 말단 염기 | 3' 말단 염기가 **G 또는 C** |
| Tm | **≥ 78 °C** |
| Tm 공식 | `Tm = 81.5 + 0.41(%GC) − 675/N − %mismatch` |

조건을 만족하는 가장 짧은 primer 쌍을 자동으로 선택합니다.

---

## E. coli 코돈 최적화

AA 변이 모드에서는 목표 아미노산에 해당하는 모든 코돈 중  
**E. coli K-12 MG1655 사용 빈도(per 1000 codons)** 가 가장 높은 코돈을 자동 선택합니다.

- 출처: Codon Usage Database (Kazusa), E. coli K-12 MG1655
- HTML 버전에서는 결과 화면의 **Codon selection** 토글을 열면 모든 후보 코돈과 빈도를 확인할 수 있습니다.

---

## 출력 예시

```
  FORWARD PRIMER   5'->3'  (sense strand)
  Primer name : MaEgtB_G100A_F
  ----------------------------------------------------------------
  Sequence : 5'-GCCTGATCGTtGCCATGCGT-3'
  Length   : 20 bp     Tm: 79.34 degC     %GC: 55.0%
  5' flank : 10 nt  |  3' flank : 8 nt
  Template : nt 289 - 308
```

---

## 주의사항

- Template 서열은 **변이 위치 양쪽으로 30 nt 이상** 확보할 것을 권장합니다.  
  여유 서열이 부족하면 조건을 만족하는 primer를 찾지 못할 수 있습니다.
- 다중 변이 Combined 설계 시 변이끼리 서열이 겹치면 설계가 불가능합니다.
- HTML 버전의 Excel 다운로드는 인터넷 연결이 필요합니다 (SheetJS CDN).
