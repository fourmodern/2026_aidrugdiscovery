"""ESMFold2 inversion 설계 — 공개 릴리스 호환 최소 어댑테이션.
근거:
 (1) disto_cond/disto_cond_mask = distogram conditioning. 저장소 테스트
     (tests/models/esmfold2_inputs_test.py FAST_PATH_OMITS)가 "inference does not use"
     라고 명시하고, prepare_input은 conditioning 미지정 시 전부 0/False 반환 → 제거는 no-op.
 (2) calculate_confidence = 공개 forward에 없는 인자. confidence는
     experimental.py forward 내부에서 `if self.confidence_head is not None:` 로
     자동 계산됨 → 인자만 제거하면 값은 그대로 나옴.
알고리즘 자체는 수정하지 않는다.
"""
import sys, time, json
sys.path.insert(0, "/workspace/esm/cookbook/tutorials")
import binder_design as bd

# --- 어댑테이션: model(**inputs, ...) 호출에서 미지원 인자만 제거 ---
_orig_fold = bd.fold_and_get_distogram
def patched_fold(model, *a, **kw):
    class _Shim:
        def __init__(self, m): self._m = m
        def __getattr__(self, n): return getattr(self._m, n)
        def __call__(self, **kwargs):
            kwargs.pop("calculate_confidence", None)   # confidence head가 자동 계산
            kwargs.pop("disto_cond", None)             # 미사용(전부 0)
            kwargs.pop("disto_cond_mask", None)        # 미사용(전부 False)
            return self._m(**kwargs)
    return _orig_fold(_Shim(model), *a, **kw)
bd.fold_and_get_distogram = patched_fold

t0=time.time()
app = bd.ESMFold2Design()
t1=time.time(); app.load(False); print("MODEL_LOAD_SEC", round(time.time()-t1,1), flush=True)
t2=time.time()
seq, traj, results = app.design(
    target_name="pd-l1", binder_name="trastuzumab_framework_vhvl",
    target_sequence=None, binder_sequence=None, is_antibody=True,
    seed=0, batch_size=1, target_hotspot_ids=None, epitope_contact_distance=12.0)
print("DESIGN_SEC", round(time.time()-t2,1), flush=True)
print("DESIGNED_SEQ:", seq, flush=True)
print("TRAJ_LEN:", len(traj), flush=True)
keys=("final_loss","iptm","ptm","distogram_iptm_proxy","cdr_distogram_iptm_proxy","critic_name","designed_sequence")
out=[{k: (float(v) if isinstance(v,(int,float)) else str(v)[:200]) for k,v in r.items() if k in keys} for r in results]
print("RESULTS_JSON:", json.dumps(out, ensure_ascii=False)[:2500], flush=True)
json.dump({"designed_sequence":seq,"traj_len":len(traj),"results":out},
          open("/workspace/esmfold2_design_results.json","w"), indent=1, ensure_ascii=False)
print("DESIGN_DONE", flush=True)
