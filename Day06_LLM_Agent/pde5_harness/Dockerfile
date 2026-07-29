FROM condaforge/mambaforge:latest

LABEL description="PDE5 inhibitor discovery agent harness (Claude Code, B-demo)"

ENV DEBIAN_FRONTEND=noninteractive
ENV BIN_DIR=/opt/bin
ENV PATH="${BIN_DIR}:${PATH}"

# System tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p ${BIN_DIR}

# Core conda packages: RDKit + requests + ChEMBL client (경량)
RUN mamba install -y -c conda-forge \
    python=3.10 \
    rdkit \
    requests \
    chembl_webresource_client \
    && mamba clean -afy

# --- 옵션(도킹 단계): smina 바이너리. 실패해도 빌드는 계속(graceful) ---
# RUN wget -q -L -O ${BIN_DIR}/smina \
#     "https://sourceforge.net/projects/smina/files/smina.static/download" \
#     && chmod +x ${BIN_DIR}/smina \
#     || echo "WARNING: smina install skipped (docking optional)"
# --- 옵션(수용체/리간드 준비): ADFRsuite prepare_receptor/prepare_ligand ---
# ADFRsuite는 https://ccsb.scripps.edu/adfr/ 에서 설치(도킹 사용 시).

WORKDIR /workspace
COPY CLAUDE.md README.md requirements.txt /workspace/
COPY scripts/ /workspace/scripts/
COPY .claude/ /workspace/.claude/

VOLUME /workspace/outputs

# 기본: 표적 조사 스크립트 실행(데모). Claude Code는 이 폴더를 열어 사용.
CMD ["python", "scripts/target_lookup.py"]
