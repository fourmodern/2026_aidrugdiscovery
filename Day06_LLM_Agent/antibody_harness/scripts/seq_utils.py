"""항체 서열 공용 유틸 — 파싱·검증·사슬 분류·가변영역(V-domain) 경계 탐지.

무-날조 원칙:
- 여기 있는 모티프 규칙은 **문헌에 공개된 보존 서열 패턴**(면역글로불린 프레임워크의
  보존 Cys 와 J-region 모티프)에 근거한 **근사 휴리스틱**이다.
- 어떤 함수도 서열을 만들어내지 않는다. 판정 불가 시 None 을 반환한다.
- 정확한 IMGT/Kabat 번호매김이 필요하면 `anarci` 또는 `abnumber` 를 설치해야 하며,
  이 모듈의 휴리스틱 결과는 `method="heuristic (approximate)"` 로만 보고해야 한다.

참고 (공개 문헌 사실, 수치가 아닌 서열 패턴):
- Ig V-domain 은 두 개의 보존 Cys 사이에 β-샌드위치를 형성한다 (Kabat H22/H92, L23/L88 근방).
- FR2 는 중쇄에서 `W-V/I-R-Q`, 경쇄에서 `W-Y/F-Q-Q` 로 시작한다.
- FR4(J-region) 는 중쇄 `W-G-x-G-T`, 경쇄 `F-G-x-G-T` 이다.
  → 이 두 모티프가 중쇄/경쇄 판별의 1차 근거가 된다.
"""
from __future__ import annotations

import re

# --- 보존 모티프 (문헌 기반 정규식) ----------------------------------------
# 중쇄 FR4 (J-region): WGxGT 계열
HEAVY_FR4 = re.compile(r"WG[QRKAEHSMGL]G[TSA]")
# 경쇄 FR4 (J-region): FGxGT 계열
LIGHT_FR4 = re.compile(r"FG[QGSTPAERK]G[TS]")
# 중쇄 FR1 끝 Cys ~ FR2 시작 Trp (그 사이가 Kabat CDR-H1 을 포함)
HEAVY_C1_W = re.compile(r"C.{8,20}?W[VIAFGLM][RKQGNHS]Q")
# 경쇄 FR1 끝 Cys ~ FR2 시작 Trp (그 사이가 Kabat CDR-L1)
LIGHT_C1_W = re.compile(r"C.{7,20}?W[YFLHVA][QLRE][QKHE]")
# FR3 시작 (중쇄): RFTIS / RVTMT / KATLT 계열
HEAVY_FR3 = re.compile(r"[RK][FVLIATM][TSVA][IFMLV][ST]")
# FR3 끝의 보존 Cys 모티프 (Kabat H92 / L88): YYC·YFC·FYC·HYC 계열
FR3_END_CYS = re.compile(r"[YFHVA][YFHCVL]C")

# --- 아미노산 스케일 (Kyte-Doolittle hydropathy, 공개 문헌 값) --------------
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def clean_seq(seq) -> str:
    """서열 문자열 정규화: 공백/개행/숫자 제거 + 대문자화."""
    if not isinstance(seq, str):
        return ""
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def parse_fasta(text: str) -> list:
    """간이 FASTA 파서. [{"id":..., "description":..., "sequence":...}, ...] 반환."""
    records, header, buf = [], None, []

    def flush():
        if header is not None:
            seq = clean_seq("".join(buf))
            if seq:
                parts = header.split(None, 1)
                records.append({
                    "id": parts[0] if parts else header,
                    "description": parts[1] if len(parts) > 1 else "",
                    "sequence": seq,
                })

    for line in (text or "").splitlines():
        line = line.strip()
        if line.startswith(">"):
            flush()
            header, buf = line[1:].strip(), []
        elif line:
            buf.append(line)
    flush()
    return records


def read_sequences(arg: str) -> list:
    """CLI 인자 하나를 서열 목록으로 해석.

    - 파일 경로(.fa/.fasta/.txt/존재하는 파일) → FASTA 파싱
    - 그 외 → 원시 서열 문자열 하나
    """
    import os

    if arg and os.path.isfile(arg):
        with open(arg, "r", encoding="utf-8") as fh:
            text = fh.read()
        if ">" in text:
            return parse_fasta(text)
        seq = clean_seq(text)
        return [{"id": os.path.basename(arg), "description": "", "sequence": seq}] if seq else []
    seq = clean_seq(arg)
    return [{"id": "input", "description": "", "sequence": seq}] if seq else []


# --- 사슬 분류 --------------------------------------------------------------
def classify_chain(seq: str) -> dict:
    """서열이 항체 중쇄/경쇄/비항체인지 보존 모티프로 판별.

    반환: {"chain_type": "heavy"|"light"|"unknown", "evidence": [...],
           "fr4_index": int|None, "method": "motif heuristic"}
    """
    s = clean_seq(seq)
    evidence, chain_type, fr4_idx = [], "unknown", None

    heavy_hits = list(HEAVY_FR4.finditer(s))
    light_hits = list(LIGHT_FR4.finditer(s))
    heavy_cw = HEAVY_C1_W.search(s)
    light_cw = LIGHT_C1_W.search(s)

    heavy_score = (1 if heavy_hits else 0) + (1 if heavy_cw else 0)
    light_score = (1 if light_hits else 0) + (1 if light_cw else 0)

    if heavy_score >= 2 and light_score >= 2:
        # 중쇄·경쇄 모티프가 모두 완전하게 존재 → scFv(또는 VH-VL 융합체) 가능성
        chain_type = "scfv"
        fr4_idx = heavy_hits[-1].start()
        evidence.append(
            f"heavy FR4 '{s[heavy_hits[-1].start():heavy_hits[-1].start() + 5]}' "
            f"@{heavy_hits[-1].start()} 와 light FR4 "
            f"'{s[light_hits[0].start():light_hits[0].start() + 5]}' @{light_hits[0].start()} 가 모두 검출"
        )
        evidence.append("단일 사슬에 VH·VL 이 함께 존재 → scFv/융합체로 판정 (CDR 분석 시 분할 필요)")
    elif heavy_score >= 2 and heavy_score >= light_score:
        chain_type = "heavy"
        fr4_idx = heavy_hits[-1].start()
        evidence.append(f"heavy FR4 motif '{s[fr4_idx:fr4_idx + 5]}' at index {fr4_idx}")
        evidence.append(f"FR1-Cys..FR2-Trp motif at index {heavy_cw.start()}")
    elif light_score >= 2:
        chain_type = "light"
        fr4_idx = light_hits[0].start()
        evidence.append(f"light FR4 motif '{s[fr4_idx:fr4_idx + 5]}' at index {fr4_idx}")
        evidence.append(f"FR1-Cys..FR2-Trp motif at index {light_cw.start()}")
    else:
        if heavy_hits:
            evidence.append("heavy FR4 motif만 존재 (FR1/FR2 모티프 불일치)")
        if light_hits:
            evidence.append("light FR4 motif만 존재 (FR1/FR2 모티프 불일치)")
        if not evidence:
            evidence.append("항체 가변영역 보존 모티프 미검출")

    return {
        "chain_type": chain_type,
        "evidence": evidence,
        "fr4_index": fr4_idx,
        "method": "motif heuristic (conserved Ig FR1-Cys/FR2-Trp/FR4 J-region patterns)",
    }


def variable_domain(seq: str) -> dict:
    """가변영역(V-domain) 대략 경계를 반환.

    start = FR1 의 보존 Cys 앞 프레임워크 시작 추정, end = FR4 모티프 + 11 (J 세그먼트 끝 근사).
    판정 불가 시 {"start": None, "end": None, "sequence": None}.
    """
    s = clean_seq(seq)
    cls = classify_chain(s)
    if cls["chain_type"] in ("unknown", "scfv") or cls["fr4_index"] is None:
        return {"start": None, "end": None, "sequence": None, "chain_type": cls["chain_type"]}

    cw = (HEAVY_C1_W if cls["chain_type"] == "heavy" else LIGHT_C1_W).search(s)
    # FR1 은 보존 Cys 앞 ~22(중쇄)/~23(경쇄) 잔기. Cys 위치에서 역산하되 음수 방지.
    offset = 22 if cls["chain_type"] == "heavy" else 23
    start = max(0, cw.start() - offset) if cw else 0
    end = min(len(s), cls["fr4_index"] + 11)
    return {
        "start": start,
        "end": end,
        "sequence": s[start:end],
        "chain_type": cls["chain_type"],
    }


# --- CDR 휴리스틱 (Kabat 근사) ----------------------------------------------
def _span(seq: str, a: int, b: int):
    if a is None or b is None or a < 0 or b > len(seq) or b <= a:
        return None
    return seq[a:b]


def _fr3_end_cys(seq: str, fr4_index: int, warnings: list):
    """CDR3 앞의 보존 Cys 인덱스를 찾는다.

    1차: 보존 모티프 [YFHVA][YFHCVL]C 의 마지막 출현(FR4 앞, 최소 3 잔기 간격).
    2차(폴백): FR4 앞 마지막 Cys — CDR3 안에 Cys 가 있으면 틀릴 수 있으므로 경고를 남긴다.
    """
    best = None
    for m in FR3_END_CYS.finditer(seq[:fr4_index]):
        c = m.end() - 1
        if fr4_index - c >= 3:
            best = c
    if best is not None:
        return best
    fallback = seq.rfind("C", 0, fr4_index)
    if fallback != -1:
        warnings.append(
            "보존 Cys 모티프([YFHVA][YFHCVL]C) 미검출 → FR4 앞 마지막 Cys 로 폴백. "
            "CDR3 내부에 Cys 가 있으면 경계가 틀릴 수 있음."
        )
        return fallback
    return None


def cdrs_heavy(seq: str) -> dict:
    """중쇄 CDR 근사 추출 (Kabat 유사 규칙, 휴리스틱).

    - CDR-H1: FR1 보존 Cys +9 잔기 ~ FR2 Trp 직전
    - CDR-H2: FR2 Trp +14 잔기 ~ FR3 시작 모티프([RK][FVLIATM][TSVA][IFMLV][ST]) 직전
    - CDR-H3: FR3 끝 보존 Cys +3 잔기 ~ FR4 (WGxGT) 직전
    실패한 CDR 은 None (날조 금지).
    """
    s = clean_seq(seq)
    out = {"CDR-H1": None, "CDR-H2": None, "CDR-H3": None, "warnings": []}

    cw = HEAVY_C1_W.search(s)
    if not cw:
        out["warnings"].append("FR1-Cys..FR2-Trp 모티프 미검출 → H1/H2 추출 불가")
    else:
        c1 = cw.start()
        w = cw.end() - 4  # 'W' 위치 (W x x Q 의 첫 글자)
        out["CDR-H1"] = _span(s, c1 + 9, w)
        h2_start = w + 14
        m3 = HEAVY_FR3.search(s, h2_start + 10, h2_start + 32)
        if m3:
            out["CDR-H2"] = _span(s, h2_start, m3.start())
        else:
            out["warnings"].append("FR3 시작 모티프 미검출 → H2 경계 불확실")

    fr4 = None
    for m in HEAVY_FR4.finditer(s):
        fr4 = m.start()
    if fr4 is None:
        out["warnings"].append("FR4(WGxGT) 미검출 → H3 추출 불가")
    else:
        c2 = _fr3_end_cys(s, fr4, out["warnings"])
        if c2 is None:
            out["warnings"].append("FR3 끝 보존 Cys 미검출 → H3 추출 불가")
        else:
            out["CDR-H3"] = _span(s, c2 + 3, fr4)
    return out


def cdrs_light(seq: str) -> dict:
    """경쇄 CDR 근사 추출 (Kabat 유사 규칙, 휴리스틱).

    - CDR-L1: FR1 보존 Cys +1 잔기 ~ FR2 Trp 직전
    - CDR-L2: FR2 Trp +15 잔기부터 7 잔기 (Kabat L50-56 고정 길이 근사)
    - CDR-L3: FR3 끝 보존 Cys +1 잔기 ~ FR4 (FGxGT) 직전
    """
    s = clean_seq(seq)
    out = {"CDR-L1": None, "CDR-L2": None, "CDR-L3": None, "warnings": []}

    cw = LIGHT_C1_W.search(s)
    if not cw:
        out["warnings"].append("FR1-Cys..FR2-Trp 모티프 미검출 → L1/L2 추출 불가")
    else:
        c1 = cw.start()
        w = cw.end() - 4
        out["CDR-L1"] = _span(s, c1 + 1, w)
        out["CDR-L2"] = _span(s, w + 15, w + 22)
        out["warnings"].append("CDR-L2 는 Kabat L50-56 고정 길이(7)를 가정한 근사값")

    m4 = LIGHT_FR4.search(s)
    if not m4:
        out["warnings"].append("FR4(FGxGT) 미검출 → L3 추출 불가")
    else:
        fr4 = m4.start()
        c2 = _fr3_end_cys(s, fr4, out["warnings"])
        if c2 is None:
            out["warnings"].append("FR3 끝 보존 Cys 미검출 → L3 추출 불가")
        else:
            out["CDR-L3"] = _span(s, c2 + 1, fr4)
    return out


def split_scfv(seq: str):
    """scFv 를 VL/VH 두 조각으로 분할. (segments, order) 반환.

    segments: [{"chain_type": "light"|"heavy", "start": int, "end": int, "sequence": str}, ...]
    분할 불가 시 (None, None).
    """
    s = clean_seq(seq)
    lm = LIGHT_FR4.search(s)
    hm = None
    for m in HEAVY_FR4.finditer(s):
        hm = m
    if not lm or not hm:
        return None, None
    l_end, h_end = lm.start(), hm.start()
    if l_end < h_end:  # VL - linker - VH
        cut = min(len(s), l_end + 11)
        return ([{"chain_type": "light", "start": 0, "end": cut, "sequence": s[:cut]},
                 {"chain_type": "heavy", "start": cut, "end": len(s), "sequence": s[cut:]}],
                "VL-linker-VH")
    cut = min(len(s), h_end + 11)  # VH - linker - VL
    return ([{"chain_type": "heavy", "start": 0, "end": cut, "sequence": s[:cut]},
             {"chain_type": "light", "start": cut, "end": len(s), "sequence": s[cut:]}],
            "VH-linker-VL")


def cdr_positions(seq: str, cdrs: dict) -> dict:
    """각 CDR 의 0-based [start, end) 인덱스를 계산(본문 위치 표기·liability 매핑용)."""
    s = clean_seq(seq)
    pos = {}
    cursor = 0
    for name in ("CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"):
        frag = cdrs.get(name)
        if not frag:
            continue
        idx = s.find(frag, cursor)
        if idx == -1:
            idx = s.find(frag)
        if idx != -1:
            pos[name] = [idx, idx + len(frag)]
            cursor = idx + len(frag)
    return pos


def try_anarci(seq: str, scheme: str = "imgt"):
    """anarci / abnumber 가 설치돼 있으면 정식 번호매김으로 CDR 을 추출.

    반환: (cdrs_dict, method_str) 또는 (None, None). 절대 값을 만들어내지 않는다.
    """
    s = clean_seq(seq)
    # 1) abnumber (ANARCI 래퍼, 사용 편의)
    try:
        from abnumber import Chain as AbChain  # type: ignore

        ch = AbChain(s, scheme=scheme)
        prefix = "H" if ch.chain_type == "H" else "L"
        cdrs = {
            f"CDR-{prefix}1": ch.cdr1_seq or None,
            f"CDR-{prefix}2": ch.cdr2_seq or None,
            f"CDR-{prefix}3": ch.cdr3_seq or None,
        }
        return cdrs, f"abnumber/ANARCI ({scheme.upper()} numbering)"
    except Exception:  # noqa: BLE001  (미설치 또는 번호매김 실패)
        pass
    # 2) anarci 직접
    try:
        from anarci import run_anarci  # type: ignore

        numbered, details, _ = run_anarci([("q", s)], scheme=scheme)
        if not numbered or not numbered[0]:
            return None, None
        dom = numbered[0][0][0]
        chain_type = details[0][0]["chain_type"]
        prefix = "H" if chain_type == "H" else "L"
        # IMGT CDR 경계 (공식 정의): CDR1 27-38, CDR2 56-65, CDR3 105-117
        bounds = {1: (27, 38), 2: (56, 65), 3: (105, 117)}
        cdrs = {}
        for k, (lo, hi) in bounds.items():
            frag = "".join(aa for (num, _ins), aa in dom if lo <= num <= hi and aa != "-")
            cdrs[f"CDR-{prefix}{k}"] = frag or None
        return cdrs, f"ANARCI ({scheme.upper()} numbering)"
    except Exception:  # noqa: BLE001
        return None, None


def numbering_backend() -> str:
    """사용 가능한 번호매김 백엔드 이름."""
    try:
        import abnumber  # noqa: F401

        return "abnumber"
    except Exception:  # noqa: BLE001
        pass
    try:
        import anarci  # noqa: F401

        return "anarci"
    except Exception:  # noqa: BLE001
        return "heuristic"
