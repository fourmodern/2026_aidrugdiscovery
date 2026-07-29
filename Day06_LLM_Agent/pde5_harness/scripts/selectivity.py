"""PDE5 vs PDE6 선택성 평가 — 정직한 한계 명시. 무-날조: 수치 예측 금지.

이 스크립트는 정량적 선택성(fold-selectivity, IC50 비)을 **예측하지 않는다**.
정량 데이터가 도구로 확보되지 않으면 "정성·확인 필요"로만 보고한다(억지 수치 금지).
제공하는 것은 (1) PDE5/PDE6 아이소자임의 기능적 맥락, (2) 교차 억제의
알려진 기전적 서술(망막 PDE6 억제 → 일시적 청녹 시각이상), (3) 실제 선택성
확인에 필요한 후속 근거(효소 어세이 / 문헌 IC50)를 명시하는 것이다.

사용: python scripts/selectivity.py [--smiles "SMILES"]
출력: 표준 결과 봉투(JSON) — provenance + verification.
"""
from __future__ import annotations
import sys, json
from verify import make_result

# PDE 아이소자임 기능 맥락(교과서적 공개 사실 — 수치 아님).
ISOFORM_CONTEXT = {
    "PDE5": {
        "gene": "PDE5A",
        "uniprot": "O76074",
        "chembl": "CHEMBL1827",
        "substrate": "cGMP",
        "tissue": "혈관 평활근, 음경 해면체, 폐혈관",
        "role": "cGMP 분해 억제 → NO-cGMP 경로 평활근 이완·혈관확장(치료 표적)",
    },
    "PDE6": {
        "gene": "PDE6A/PDE6B/PDE6C 등",
        "substrate": "cGMP",
        "tissue": "망막 광수용체(rod/cone) 외절",
        "role": "광전달(phototransduction) 캐스케이드의 효과기 효소",
    },
}

# 교차 억제 기전 서술(정성). 특정 fold-selectivity 수치는 여기에 기재하지 않음.
CROSS_INHIBITION_MECHANISM = (
    "PDE5 저해제는 PDE5와 구조적으로 유사한 망막 PDE6를 부분적으로 교차 억제할 수 있다. "
    "PDE6가 억제되면 광수용체 내 cGMP 항상성이 교란되어 일시적 색각 이상(청녹색 색조 변화, "
    "cyanopsia)·광민감성이 유발될 수 있다. 이는 PDE5 저해제 계열에서 관찰되는 대표적 "
    "표적-외(off-target) 시각 부작용의 기전적 근거로 알려져 있다."
)

# 실제 선택성 확인에 필요한 근거(무-날조: 아래를 확보하기 전에는 수치 단정 금지).
CONFIRMATION_REQUIRED = [
    "PDE5 및 PDE6 효소에 대한 동일 조건 IC50/Ki 병렬 어세이 (fold-selectivity 산출의 유일한 근거)",
    "ChEMBL 등에서 화합물별 PDE5/PDE6 활성값 병렬 조회 (target: PDE6 여러 아이소자임)",
    "문헌 보고 선택성 비의 출처(PMID/DOI) 명시 — 인용 없이 수치 사용 금지",
]

REFERENCES = [
    "Boolell M et al. Sildenafil: an orally active type 5 cyclic GMP-specific "
    "phosphodiesterase inhibitor for the treatment of penile erectile dysfunction. "
    "Int J Impot Res. 1996;8(2):47-52.",
    "Ghofrani HA, Osterloh IH, Grimminger F. Sildenafil: from angina to erectile "
    "dysfunction to pulmonary hypertension and beyond. Nat Rev Drug Discov. "
    "2006;5(8):689-702. doi:10.1038/nrd2030",
]


def main():
    smiles = None
    if "--smiles" in sys.argv:
        i = sys.argv.index("--smiles")
        if i + 1 < len(sys.argv):
            smiles = sys.argv[i + 1]

    result = {
        "quantitative_selectivity_predicted": False,
        "assessment": "정성(qualitative) — 정량 fold-selectivity 미산출",
        "isoform_context": ISOFORM_CONTEXT,
        "cross_inhibition_mechanism": CROSS_INHIBITION_MECHANISM,
        "confirmation_required": CONFIRMATION_REQUIRED,
        "references": REFERENCES,
        "input_smiles": smiles,
        "caveat": (
            "이 도구는 PDE5/PDE6 선택성을 정량 예측하지 않는다. 화합물별 선택성 결론은 "
            "위 confirmation_required 근거(병렬 IC50 어세이 또는 인용된 문헌 수치)가 확보된 "
            "뒤에만 내릴 수 있다. 수치·비율 날조 금지."
        ),
    }

    # 검증: 정량 수치를 만들어내지 않았음(정직성) + 기전/근거 서술 존재.
    checks = [
        ("정량 fold-selectivity 미날조(정성 보고)",
         result["quantitative_selectivity_predicted"] is False),
        ("교차 억제 기전(PDE6 망막) 서술 존재",
         "PDE6" in result["cross_inhibition_mechanism"]
         and "망막" in result["cross_inhibition_mechanism"]),
        ("후속 확인 근거 목록 존재", len(result["confirmation_required"]) > 0),
        ("참고문헌 실재(Boolell 1996 / Ghofrani 2006)", len(result["references"]) >= 2),
    ]
    env = make_result(
        result, "PDE5/PDE6 isoform context (public textbook facts) + literature",
        f"selectivity qualitative; smiles={smiles}", checks,
        notes=("PDE5 vs PDE6 선택성은 정성적으로만 보고. 정량 예측 없음 → 억지 수치 금지. "
               "교차 억제(망막 PDE6) → 일시적 청녹 시각이상 기전 서술. 실제 선택성은 어세이/문헌 확인 권장."),
    )
    print(json.dumps(env, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
