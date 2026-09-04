#!/usr/bin/env python3
"""학회 구두발표용 슬라이드 (16:9, 12분 · 12장).

`make_docs.py` 가 만드는 pptx 는 보고서를 기계적으로 옮긴 것이라 발표에 쓸 수 없다.
발표는 문서 낭독이 아니라 **하나의 논증**이다. 그래서 별도로 만든다.

원칙
  - 한 장에 메시지 하나. 슬라이드 제목이 그 메시지 자체다 (주제어가 아니라 문장).
  - 숫자는 크게. 뒷줄에서 읽혀야 한다.
  - 그림은 꽉 채운다. 캡션은 슬라이드 제목이 대신한다.
  - 발표자 노트에 말할 내용을 적는다 — 슬라이드에 적지 않는다.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
from pptx import Presentation
from pptx.util import Mm, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

ROOT = Path(__file__).resolve().parent.parent
SR = ROOT / "sample_run"
FIG = SR / "report" / "figures_controlled"
FONT = "Pretendard"
W, H = 338.7, 190.5                    # 16:9 (13.333 × 7.5 in)

INK    = RGBColor(0x1A, 0x1A, 0x1A)
NAVY   = RGBColor(0x0B, 0x2A, 0x4A)
BLUE   = RGBColor(0x00, 0x72, 0xB2)
ORANGE = RGBColor(0xD5, 0x5E, 0x00)
GREEN  = RGBColor(0x00, 0x9E, 0x73)
GREY   = RGBColor(0x5A, 0x63, 0x6B)
LGREY  = RGBColor(0xF3, 0xF5, 0xF7)
WHITE  = RGBColor(0xFF, 0xFF, 0xFF)
PALE   = RGBColor(0x9F, 0xC5, 0xE0)


def load(n):
    p = SR / n
    if not p.exists(): return None
    d = json.loads(p.read_text()); return d.get("result", d)


class Deck:
    def __init__(self):
        self.prs = Presentation()
        self.prs.slide_width, self.prs.slide_height = Mm(W), Mm(H)
        self.blank = self.prs.slide_layouts[6]
        self.n = 0

    def box(self, s, x, y, w, h, fill=None, line=None, lw=1.0,
            shape=MSO_SHAPE.ROUNDED_RECTANGLE):
        sh = s.shapes.add_shape(shape, Mm(x), Mm(y), Mm(w), Mm(h))
        if fill: sh.fill.solid(); sh.fill.fore_color.rgb = fill
        else: sh.fill.background()
        if line: sh.line.color.rgb = line; sh.line.width = Pt(lw)
        else: sh.line.fill.background()
        sh.shadow.inherit = False
        return sh

    def text(self, s, x, y, w, h, runs, align=PP_ALIGN.LEFT,
             anchor=MSO_ANCHOR.TOP, space=0):
        tb = s.shapes.add_textbox(Mm(x), Mm(y), Mm(w), Mm(h))
        tf = tb.text_frame; tf.word_wrap = True; tf.vertical_anchor = anchor
        tf.margin_left = tf.margin_right = Mm(1); tf.margin_top = tf.margin_bottom = 0
        paras = runs if isinstance(runs[0], list) else [runs]
        for i, para in enumerate(paras):
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            p.alignment = align; p.space_after = Pt(space)
            for t, sz, colr, *rest in para:
                r = p.add_run(); r.text = t
                r.font.size = Pt(sz); r.font.color.rgb = colr
                r.font.bold = bool(rest and rest[0]); r.font.name = FONT
        return tb

    def slide(self, title, kicker="", note=""):
        """제목은 주제어가 아니라 그 장에서 하려는 말 자체다."""
        s = self.prs.slides.add_slide(self.blank); self.n += 1
        self.box(s, 0, 0, W, 2.2, fill=BLUE, shape=MSO_SHAPE.RECTANGLE)
        if kicker:
            self.text(s, 14, 6, 120, 6, [(kicker, 12, BLUE, True)])
        self.text(s, 14, 12, W - 28, 16, [(title, 25, NAVY, True)])
        self.text(s, W - 26, H - 11, 16, 6, [(str(self.n), 12, GREY)],
                  align=PP_ALIGN.RIGHT)
        if note:
            s.notes_slide.notes_text_frame.text = note
        return s

    def pic(self, s, name, x, y, w, maxh=None):
        f = FIG / name
        if not f.exists():
            self.box(s, x, y, w, maxh or 60, fill=LGREY, line=GREY)
            return
        p = s.shapes.add_picture(str(f), Mm(x), Mm(y), width=Mm(w))
        if maxh and p.height / 36000 > maxh:
            r = maxh / (p.height / 36000)
            p.height = int(p.height * r); p.width = int(p.width * r)
            p.left = Mm(x + (w - p.width / 36000) / 2)
        return p

    def bignum(self, s, x, y, w, val, lab, colr=ORANGE, sz=54):
        self.text(s, x, y, w, 22, [(val, sz, colr, True)], align=PP_ALIGN.CENTER)
        self.text(s, x, y + 20, w, 10, [(lab, 13, GREY)], align=PP_ALIGN.CENTER)

    def save(self, out):
        self.prs.save(str(out)); return self.n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(SR / "report" / "docs" / "talk_16x9.pptx"))
    a = ap.parse_args()
    ds, dk, an = load("dataset_controlled.json"), load("docking_controlled.json"), \
                 load("analysis_controlled.json")
    tc, sw = load("terms_controlled.json"), load("exhaustiveness_sweep.json")
    prot, col = load("protonation_test.json"), load("collapse_diagnosis.json")
    oldp = load("old_set_new_protocol.json")
    if not (ds and dk and an):
        raise SystemExit("산출 파일이 없다 — 슬라이드를 만들지 않는다.")
    T = an["arms"]["top_pose"]; PV = T["pose_validity"]; wb = T["within_similarity_bin"]
    ACI = T["auc_ci95"]; NSEN = T["near_sensitivity"]; ctrl = dk["control_redock"]
    d = Deck()

    # 1 표지
    s = d.prs.slides.add_slide(d.blank); d.n += 1
    d.box(s, 0, 0, W, H, fill=NAVY, shape=MSO_SHAPE.RECTANGLE)
    d.text(s, 24, 52, W - 48, 30,
           [("골격 유사성과 역가를 직교화한\nPDE5A 도킹 재평가", 40, WHITE, True)],
           align=PP_ALIGN.CENTER)
    d.text(s, 24, 100, W - 48, 12,
           [("선행 상관의 재현 실패, 그리고 자세 타당도가 정하는 해석의 상한", 18, PALE)],
           align=PP_ALIGN.CENTER)
    d.text(s, 24, 132, W - 48, 20,
           [[(f"smina / AutoDock Vina 1.1.2  ·  PDB 1UDT  ·  ChEMBL PDE5A  n = {an['n']}",
              13, PALE)],
            [("전 산출물·코드·리뷰 기록 공개 — github.com/fourmodern/2026_aidrugdiscovery",
              13, WHITE, True)]], align=PP_ALIGN.CENTER, space=4)
    s.notes_slide.notes_text_frame.text = (
        "12분 발표. 결론을 3번 슬라이드에서 먼저 말하고, 나머지는 그 근거를 쌓는 구성입니다.")

    # 2 문제
    s = d.slide("도킹 벤치마크는 두 가지를 동시에 재고 있다", "문제",
                "표적 구조는 대개 holo 입니다. 그 리간드를 닮은 화합물은 점수도 좋고 실제로도 "
                "강력합니다. 두 경로가 같은 방향이라 상관은 나오는데, 그게 채점 능력인지 "
                "골격 닮음인지 구분이 안 됩니다.")
    d.pic(s, "fig01_graphical_abstract.png", 12, 36, W - 24, maxh=132)

    # 3 결론 먼저
    s = d.slide("결론부터 — 골격을 분리하니 상관이 사라졌다", "핵심",
                "여기서 결론을 먼저 말합니다. 뒤는 전부 이 숫자를 어떻게 얻었고 무엇을 "
                "배제했는지입니다.")
    d.box(s, 20, 46, W - 40, 46, fill=RGBColor(0xFD, 0xF0, 0xE6), line=ORANGE, lw=2)
    d.text(s, 26, 55, W - 52, 14,
           [(f"ρ = {T['spearman']:+.3f}    (95% CI {T['ci95'][0]:+.2f} ~ "
             f"{T['ci95'][1]:+.2f},   p = {T['perm_p']})", 34, ORANGE, True)],
           align=PP_ALIGN.CENTER)
    d.text(s, 26, 74, W - 52, 12,
           [(f"n = 30 에서 얻었던 ρ = {col['old_design']['spearman_top_pose']:+.3f} 는 "
             f"재현되지 않았다", 17, INK)], align=PP_ALIGN.CENTER)
    for i, (v, lab, c) in enumerate([
            (f"{tc['models']['single_vina_score']['Q2_loo']:+.3f}" if tc else "—",
             "교차검증 Q²  (기준 0.3)", ORANGE),
            (f"{ACI['all']['auc']}", "ROC-AUC 전체  (기준 0.7)", ORANGE),
            (f"{1 - PV['all']['frac_under_threshold']:.0%}", "자세가 2 Å 기준 밖", ORANGE)]):
        d.bignum(s, 20 + i * ((W - 40) / 3), 112, (W - 40) / 3, v, lab, c, sz=50)

    # 4 설계
    s = d.slide("역가와 골격 유사도를 9칸 격자로 직교화했다", "설계",
                "ChEMBL 전수를 받아 화합물당 중앙값으로 집계하고, 역가 3구간과 실데나필 "
                "Tanimoto 3구간의 격자에서 균등 추출했습니다. 교란 상관이 0.281에서 "
                "0.175로 떨어집니다.")
    d.pic(s, "fig02_dataset.png", 12, 40, W - 24, maxh=88)
    d.text(s, 20, 136, W - 40, 18,
           [[(f"ChEMBL PDE5A 활성 {ds['records_scanned']} 레코드 전수 → 화합물당 중앙값 "
              f"({ds['pool_size']}종) → 9칸 격자 균등 추출 n = {an['n']}", 15, INK)],
            [(f"교란 상관 Pearson "
              f"{ds['confound_pearson_r_potency_vs_similarity']:+.3f}  "
              f"(역가로만 층화한 옛 설계 +0.281)", 16, GREEN, True)]],
           align=PP_ALIGN.CENTER, space=3)

    # 5 대조
    s = d.slide("재도킹 대조는 하나가 아니라 둘로 쪼개야 한다", "방법",
                "C1은 상위 모드 중 최선, C2는 점수 1위 자세입니다. 합쳐서 하나로 보고하면 "
                "탐색 실패인지 채점 실패인지 알 수 없습니다.")
    d.pic(s, "fig03_control.png", 12, 38, W - 24, maxh=94)
    for i, (nm, v, ok_, c) in enumerate([
            ("C1  샘플링", f"{ctrl['rmsd_best_of_modes_angstrom']:.2f} Å", "PASS", GREEN),
            ("C2  채점", f"{ctrl['rmsd_top_pose_angstrom']:.2f} Å", "FAIL", ORANGE)]):
        x = 60 + i * 110
        d.box(s, x, 138, 100, 22, fill=c, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
        d.text(s, x + 6, 142, 88, 14,
               [(f"{nm}   {v}   {ok_}", 19, WHITE, True)], align=PP_ALIGN.CENTER)

    # 6 스윕
    s = d.slide("탐색을 32배 늘려도 채점은 미동도 없었다", "원인 배제 ①",
                "깊이 4에서 128까지. C1은 통과율 0에서 100%로 고쳐지는데 C2는 전 구간 0%이고 "
                "최고 점수는 아예 변하지 않습니다. 탐색은 얕은 깊이에서 이미 전역 최소를 "
                "찾았고 그 최소가 틀린 자리에 있습니다.")
    d.pic(s, "fig08_exhaustiveness.png", 12, 38, W - 24, maxh=104)
    d.text(s, 20, 148, W - 40, 12,
           [[("C1 통과율 0% → 100%", 17, GREEN, True), ("      ", 17, INK),
             (f"C2 전 구간 0%, 변동 {sw['diagnosis']['c2_range_angstrom']} Å", 17, ORANGE, True),
             ("      ", 17, INK),
             (f"최고 점수 변동 {sw['diagnosis']['score_range_kcal']} kcal/mol", 17, ORANGE, True)]],
           align=PP_ALIGN.CENTER)

    # 7 프로토네이션
    s = d.slide("프로토네이션도 원인이 아니었다 — 그리고 그럴 수 없었다", "원인 배제 ②",
                "수용체는 pH 7.4인데 리간드는 중성이라는 비대칭이 있었습니다. 맞춰도 점수가 "
                "그대로입니다. 이유는 기본 Vina 함수에 정전기 항이 없기 때문입니다.")
    rows = [("수용체 pH7.4 · 리간드 pH7.4  (옳은 짝)", "recpH74_ligpH74"),
            ("수용체 pH7.4 · 리간드 중성  (현행)", "recpH74_ligneutral"),
            ("수용체 중성 · 리간드 중성", "recneutral_ligneutral"),
            ("수용체 중성 · 리간드 pH7.4", "recneutral_ligpH74")]
    by = {x["condition"]: x for x in (prot or {}).get("summary", [])}
    d.text(s, 22, 40, 150, 8, [("조건", 14, GREY, True)])
    for j, lab in enumerate(("C1 (Å)", "C2 (Å)", "점수", "C2 통과")):
        d.text(s, 180 + j * 38, 40, 36, 8, [(lab, 14, GREY, True)], align=PP_ALIGN.CENTER)
    for i, (lab, key) in enumerate(rows):
        r = by.get(key, {}); yy = 50 + i * 15
        d.box(s, 20, yy, W - 40, 13.5, fill=(LGREY if i % 2 == 0 else WHITE),
              shape=MSO_SHAPE.RECTANGLE)
        d.text(s, 24, yy + 1.5, 154, 11, [(lab, 15, INK, i == 0)],
               anchor=MSO_ANCHOR.MIDDLE)
        for j, v in enumerate((f"{r.get('c1_best_rmsd_mean', 0):.2f}",
                               f"{r.get('c2_top_rmsd_mean', 0):.2f}",
                               f"{r.get('top_score_mean', 0):.2f}",
                               f"{r.get('c2_pass_rate', 0):.0%}")):
            d.text(s, 180 + j * 38, yy + 1.5, 36, 11,
                   [(v, 16, ORANGE if j == 3 else INK, True)],
                   align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
    d.box(s, 20, 118, W - 40, 26, fill=RGBColor(0xFD, 0xF0, 0xE6), line=ORANGE, lw=1.4)
    d.text(s, 26, 123, W - 52, 18,
           [[("기본 Vina 채점 함수에는 정전기 항이 없다.", 18, ORANGE, True)],
            [("형식전하는 점수에 닿지 못하고, 바뀐 수소결합 타이핑도 그 질소가 결합 거리 "
              "밖이라 항 값이 움직이지 않았다.", 15, INK)]], space=3)

    # 8 붕괴 원인
    s = d.slide("무너진 원인은 골격도, 설계도, 프로토콜도 아니었다", "원인 배제 ③",
                "초안은 골격 교란 때문이라고 적었는데, 옛 데이터에 같은 통제를 적용하니 "
                "1.2%만 줄었습니다. 설계 재현도, 프로토콜 재도킹도 원인이 아니었습니다. "
                "남은 것은 화합물 집합 자체입니다.")
    tests = [("① 옛 데이터에 골격 통제",
              f"{col['old_design']['spearman_top_pose']:+.3f} → "
              f"{col['old_design']['partial_controlling_tanimoto']:+.3f}",
              f"감쇠 {col['old_design']['attenuation_from_scaffold_control']:.1%}"),
             ("② 옛 데이터의 점수 vs Tanimoto",
              f"{col['old_design']['spearman_score_vs_tanimoto']:+.3f}",
              f"p = {col['old_design']['perm_p_score_vs_tanimoto']}"),
             ("③ 옛 설계를 새 데이터에 재현",
              f"중앙값 {col['test_1_design']['median']:+.3f}",
              f"−0.4 미만 {col['test_1_design']['frac_below_-0.4']:.1%}"),
             ("④ 옛 30건을 새 프로토콜로 재도킹",
              f"{oldp['rho_old_protocol']:+.3f} → {oldp['rho_new_protocol']:+.3f}" if oldp else "—",
              f"점수 상관 {oldp['spearman_between_protocols']:+.3f}" if oldp else "—")]
    for i, (nm, v1, v2) in enumerate(tests):
        yy = 42 + i * 20
        d.box(s, 20, yy, W - 40, 18, fill=(LGREY if i % 2 == 0 else WHITE),
              shape=MSO_SHAPE.RECTANGLE)
        d.text(s, 25, yy + 2, 130, 14, [(nm, 16, INK)], anchor=MSO_ANCHOR.MIDDLE)
        d.text(s, 158, yy + 2, 82, 14, [(v1, 18, ORANGE, True)],
               align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
        d.text(s, 244, yy + 2, 74, 14, [(v2, 15, GREY)],
               align=PP_ALIGN.RIGHT, anchor=MSO_ANCHOR.MIDDLE)
    d.text(s, 20, 128, W - 40, 16,
           [("네 검정이 모두 후보를 배제한다. 원인은 그 화합물 집합 자체이며, "
             "어떤 성질 때문인지는 이 데이터로 더 좁히지 못했다.", 17, NAVY, True)],
           align=PP_ALIGN.CENTER)

    # 9 통제 후 상관
    s = d.slide("교란을 통제해도, 대역을 나눠도 신호는 없다", "결과 ①",
                "골격과 크기를 각각·동시에 통제했습니다. far와 mid 대역은 0 근처이고 "
                "near만 유의합니다.")
    d.pic(s, "fig05_forest.png", 60, 34, W - 120, maxh=118)
    d.text(s, 20, 156, W - 40, 12,
           [[("전체 ", 16, GREY), (f"{T['spearman']:+.3f}", 18, ORANGE, True),
             ("    골격 통제 ", 16, GREY),
             (f"{T['partial_spearman_controlling_tanimoto']:+.3f}", 18, ORANGE, True),
             ("    골격+크기 ", 16, GREY),
             (f"{T['partial_spearman_controlling_both']:+.3f}", 18, ORANGE, True)]],
           align=PP_ALIGN.CENTER)

    # 10 near 대역
    s = d.slide("한 대역만 유의하지만, 그 신호는 칸 하나에 얹혀 있다", "결과 ②",
                "near 대역만 유의하고 보정 후에도 살아남습니다. 그런데 그 대역의 약함 칸이 "
                "전수라 민감도 분석이 안 되고, AUC 신뢰구간이 기준 0.7을 포함합니다. "
                "가설로만 남깁니다.")
    for i, (v, lab, c) in enumerate([
            (f"{wb['near']['spearman']:+.3f}", f"near 대역 ρ  (p = {wb['near']['perm_p']})", GREEN),
            (f"{ACI['near']['auc']}",
             f"AUC  95% CI [{ACI['near']['ci95'][0]:.2f}, {ACI['near']['ci95'][1]:.2f}]", ORANGE),
            (f"{NSEN['drop_weak_cell']['spearman']:+.3f}", "전수 칸 제외 시 ρ", ORANGE)]):
        d.bignum(s, 16 + i * ((W - 32) / 3), 56, (W - 32) / 3, v, lab, c, sz=48)
    d.box(s, 24, 108, W - 48, 40, fill=RGBColor(0xFD, 0xF0, 0xE6), line=ORANGE, lw=1.6)
    d.text(s, 30, 114, W - 60, 30,
           [[("이 결과에 무게를 실을 수 없다", 20, ORANGE, True)],
            [("(a) AUC 신뢰구간이 사전 기준 0.7 을 포함한다   "
              "(b) 전수(census)인 약함 칸을 빼면 유의하지 않다   "
              "(c) 그 칸의 역가가 어세이 바닥값에 몰려 있다", 15, INK)]], space=4)

    # 11 자세 타당도
    s = d.slide("우리가 점수를 매긴 자세의 81%는 결정 자세 틀 밖에 있다", "해석의 상한",
                "이게 이 연구에서 가장 큰 제약입니다. 신호 없음이라는 진술이 대부분 틀린 "
                "자세에 매긴 점수에 대한 것이라는 뜻입니다.")
    for i, (v, lab) in enumerate([
            (f"{PV['all']['median_rmsd']:.2f} Å", "1위 자세 MCS-RMSD 중앙값"),
            (f"{PV['all']['frac_under_threshold']:.1%}", "2 Å 기준 통과 (전체)"),
            (f"{PV['far']['frac_under_threshold']:.1%}", "far 대역 통과율"),
            (f"{PV['near']['frac_under_threshold']:.1%}", "near 대역 통과율")]):
        d.bignum(s, 12 + i * ((W - 24) / 4), 56, (W - 24) / 4, v, lab, ORANGE, sz=42)
    d.text(s, 24, 112, W - 48, 40,
           [[("이 지표가 무엇을 재는지 분명히 해 둔다.", 17, NAVY, True)],
            [("실데나필 외의 화합물에는 실험 결정 자세가 없다. 이 값은 “자세가 맞는가”가 "
              "아니라 “공유 부분구조가 공결정 리간드 위에 겹치는가”를 잰다. 골격이 다른 "
              "far 대역에서는 해석이 특히 약하므로 대역별로 읽어야 한다.", 15, INK)]], space=4)

    # 12 결론
    s = d.slide("정리 — 무엇을 말할 수 있고 무엇을 말할 수 없는가", "결론",
                "실무 지침 셋으로 마무리합니다.")
    items = [("정량 예측은 실패했다", "Q² 세 모델 모두 0 근처, 기준 0.3 에 미달", ORANGE),
             ("선별도 전체적으로는 작동하지 않았다", "한 대역의 신호는 칸 하나에 의존", ORANGE),
             ("네 후보를 실험으로 제거했다",
              "탐색 깊이 · 프로토네이션 · 프로토콜 · 골격 교란", GREEN),
             ("남은 후보", "물 제거 · 강체 수용체 · 단일 배좌 · holo 편향 · 채점 함수", GREY)]
    for i, (h_, sub, c) in enumerate(items):
        yy = 40 + i * 21
        d.box(s, 20, yy, 4, 17, fill=c, shape=MSO_SHAPE.RECTANGLE)
        d.text(s, 29, yy, W - 50, 10, [(h_, 19, c if c != GREY else NAVY, True)])
        d.text(s, 29, yy + 9.5, W - 50, 8, [(sub, 14, GREY)])
    d.box(s, 20, 128, W - 40, 32, fill=LGREY, line=NAVY, lw=1.6)
    d.text(s, 26, 133, W - 52, 24,
           [[("실무 지침", 17, NAVY, True)],
            [("① 재도킹 대조를 C1/C2 로 쪼개고 깊이를 흔들어 병목을 특정하라   "
              "② 점수를 매긴 자세가 어디에 놓였는지 먼저 보고하라   "
              "③ 벤치마크의 골격 분포를 통제하라", 15, INK)]], space=3)

    # 13 마무리
    s = d.prs.slides.add_slide(d.blank); d.n += 1
    d.box(s, 0, 0, W, H, fill=NAVY, shape=MSO_SHAPE.RECTANGLE)
    d.text(s, 24, 56, W - 48, 24,
           [("같은 데이터로 세 판본, 세 결론.\n매번 게이트는 전부 통과했다.", 32, WHITE, True)],
           align=PP_ALIGN.CENTER)
    d.text(s, 24, 100, W - 48, 14,
           [("오류를 찾은 것은 자동 검증이 아니라 외부 비평이었다. 게이트는 수치가 "
             "산출되었는지를 검사할 뿐, 결론이 데이터에서 따라 나오는지는 검사하지 못한다.",
             16, PALE)], align=PP_ALIGN.CENTER)
    d.text(s, 24, 138, W - 48, 16,
           [[("전 산출물 · 코드 · 리뷰 기록 공개", 15, PALE)],
            [("github.com/fourmodern/2026_aidrugdiscovery → Day06_LLM_Agent/pde5_harness",
              17, WHITE, True)]], align=PP_ALIGN.CENTER, space=3)
    s.notes_slide.notes_text_frame.text = (
        "질문 대비: (1) 왜 PDE5A 인가 — 결정 구조와 활성 데이터가 모두 풍부해서. "
        "(2) 다른 표적에서도 같은가 — 검정하지 않았고 그것이 가장 큰 한계. "
        "(3) 물을 넣으면 달라지나 — 배제하지 못한 후보이며 다음 실험.")

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    n = d.save(out)
    print(f"{out}  ({n}장, 16:9)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
