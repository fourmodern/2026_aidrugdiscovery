#!/usr/bin/env python3
"""음성 대조 — 실제 도구와 실제 파이프라인이 실패를 막는지 관측한다.

[개정 2026-09-04] 이전 판은 손으로 만든 딕셔너리를 gate() 에 넣고 3줄짜리 모의 루프를
돌렸다. 그것은 "gate() 가 False 를 반환한다"와 "break 가 루프를 끊는다"만 보여줄 뿐,
실제 파이프라인이 멈추는지는 검증하지 못했다. 이 판은 진짜 도구를 실패시킨다.

세 가지를 관측한다.
  A. mol_properties.py 에 파싱 불가 SMILES 를 주면 실제로 passed=false 봉투가 나오는가
  B. selectivity.py 에 쓰레기 SMILES 를 주면 passed=false 가 나오는가 (수정된 술어)
  C. 실패 봉투를 파이프라인 2번 위치에 넣으면 3번이 실행되지 않는가
     (주입하는 봉투는 B 에서 selectivity.py 가 실제로 만든 것이다)
"""
import json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PY = sys.executable
results = []


def run(script, args=None, stdin_text=None):
    p = subprocess.run([PY, str(HERE / script)] + (args or []), input=stdin_text,
                       capture_output=True, text=True, cwd=ROOT)
    try:
        return json.loads(p.stdout), p.returncode
    except json.JSONDecodeError:
        return None, p.returncode


print("== A. mol_properties.py 에 파싱 불가 SMILES 주입 ==")
bad = json.dumps({"result": [{"canonical_smiles": "NOT_A_MOLECULE_@@@"},
                             {"canonical_smiles": "###"}]})
env, _ = run("mol_properties.py", ["--stdin"], bad)
a_ok = env is not None and env["verification"]["passed"] is False
print(f"  passed = {env['verification']['passed'] if env else 'None'}"
      f"  실패항목 = {[c['check'] for c in env['verification']['checks'] if not c['passed']] if env else []}")
print(f"  → {'PASS — 실제 도구가 실패를 보고함' if a_ok else 'FAIL — 실패가 보고되지 않음'}")
results.append(("A mol_properties 실패 보고", a_ok))

print("\n== B. selectivity.py 에 쓰레기 SMILES 주입 (수정된 입력 의존 술어) ==")
env, _ = run("selectivity.py", ["--smiles", "NOT_A_MOLECULE_@@@"])
b_ok = env is not None and env["verification"]["passed"] is False
print(f"  passed = {env['verification']['passed'] if env else 'None'}"
      f"  실패항목 = {[c['check'][:26] for c in env['verification']['checks'] if not c['passed']] if env else []}")
print(f"  → {'PASS — 입력 의존 술어가 발동함' if b_ok else 'FAIL — 여전히 반증 불가'}")
results.append(("B selectivity 반증 가능", b_ok))

print("\n== C. run_harness 실제 경로가 실패 봉투를 만나면 멈추는가 ==")
print("   (B 에서 얻은 실제 실패 봉투를 2번 위치에 주입하고 3번이 실행되는지 본다)")
sys.path.insert(0, str(HERE))
from verify import gate  # noqa: E402

executed = []
steps = [("target-lookup", {"verification": {"passed": True, "checks": []}}),
         ("selectivity-check(주입)", env if env else {"verification": {"passed": False, "checks": []}}),
         ("mol-properties", {"verification": {"passed": True, "checks": []}})]
for name, e in steps:
    if not gate(e, name):
        print(f"   {name} 에서 중단. 이후 단계 미실행.")
        break
    executed.append(name)
c_ok = "mol-properties" not in executed
print(f"   실행된 단계: {executed}")
print(f"  → {'PASS — 실패 봉투가 다음 단계를 막음' if c_ok else 'FAIL — 막지 못함'}")
results.append(("C 실패 봉투가 다음 단계 차단", c_ok))

print("\n" + "=" * 60)
for name, ok in results:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
allok = all(ok for _, ok in results)
print(f"\n판정: {'PASS — 실제 도구 수준에서 차단이 관측됨' if allok else 'FAIL'}")
print("한계: C 는 run_harness 의 gate 호출 규약을 재현한 것이며,")
print("      run_harness.py 프로세스 자체를 실패 입력으로 끝까지 돌린 것은 아니다.")
sys.exit(0 if allok else 1)
