#!/usr/bin/env python3
"""원본 batch_run.build (patch_build.py 적용본) 으로 <poly>+PFOA 스크린 계 하나를 만든다.
사용: (~/mfpvector 에서) python3 build_one.py PTFE
유일한 변경: packmol 이 '완전 배치 아님'(비정상 종료 코드)으로 끝나도 출력 파일이 있으면 계속 진행
(원본 실행 때 쓰인 packmol 버전은 이 경우에도 0 을 돌려줬을 수 있음; 겹침은 이어지는 진공 EM 이 푼다)."""
import sys, os, subprocess
sys.path.insert(0, os.path.expanduser("~/mfpvector/code"))
import batch_run as b

_run = b.run


def run(cmd, inp=None, shell=False, quiet=True):
    if shell and isinstance(cmd, str) and cmd.startswith("packmol"):
        r = subprocess.run(cmd, input=inp, capture_output=True, text=True, shell=True)
        inpf = cmd.split("<")[1].strip()
        out = [l.split()[1] for l in open(inpf) if l.startswith("output")][0]
        if r.returncode != 0:
            if os.path.exists(out) and os.path.getsize(out) > 0:
                b.log(f"packmol 완전 배치 아님 (code {r.returncode}) — 출력으로 계속")
                return r
            raise RuntimeError("packmol 실패:\n" + r.stdout[-1500:])
        return r
    return _run(cmd, inp=inp, shell=shell, quiet=quiet)


b.run = run
p = sys.argv[1]
print(b.build(p, "PFOA", f"systems/{p}_PFOA"))
