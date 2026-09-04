#!/usr/bin/env python3
"""실행 환경을 기계적으로 캡처해 sample_run/environment.json 에 쓴다.

보고서가 "환경 문자열도 기계 캡처값" 이라고 주장하려면 그 주장을 뒷받침하는 코드가
저장소에 있어야 한다. 리뷰에서 "그 스크립트가 없다" 는 지적을 받아 추가했다.
사람이 값을 타이핑할 여지를 남기지 않는다.
"""
import importlib.metadata as md
import json
import sys
from pathlib import Path

PACKAGES = ("rdkit", "chembl-webresource-client", "requests", "matplotlib", "numpy")


def capture() -> dict:
    env = {}
    for name in PACKAGES:
        try:
            env[name] = md.version(name)
        except md.PackageNotFoundError:
            env[name] = "미설치"
    env["python"] = sys.version.split()[0]
    return env


def main() -> int:
    out = Path(__file__).resolve().parent.parent / "sample_run" / "environment.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    env = capture()
    out.write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"{out} 갱신")
    for k, v in env.items():
        print(f"  {k:28s} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
