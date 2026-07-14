#!/usr/bin/env python
"""
main.py
~~~~~~~
Demonstration script showing the features and extensible architecture of the 
Universal CRISPR Guide RNA Design Platform - CRISPR-UNI.
Powered by Polymorphic Deep Learning Predictors.
"""

import logging

try:
    import readline
except ImportError:
    pass

from crispr_designer import (
    CRISPRDesigner, 
    CRISPRSystem, 
    CRISPR_REGISTRY, 
    ORGANISM_REGISTRY,
    PredictorFactory,
    BasePredictor
)
from crispr_designer.predictors import ModelWeightsMissingError

# Multi-language localization dictionary
LOCALIZATION = {
    "en": {
        "title_bar": "CRISPR-UNI Interactive Design Tool",
        "welcome": "Welcome to the extensible biological guide RNA design platform.",
        "exit_info": "Type 'exit' at any prompt to quit the application.",
        "step1_title": "Step 1: Enter Target Sequence",
        "step1_desc1": "  - Please enter the DNA or RNA target sequence to design guides.",
        "step1_desc2": "  - Example DNA: ATGCGATCGATCGATCGATCGATCGATCGAATCGATCGATCGATCGGGCATCGATCG",
        "step1_desc3": "  - Example RNA: AUCGAUCGAUCGAUCGAUCGAUCGAUGGC",
        "step1_prompt": "  >> Enter target sequence: ",
        "invalid_seq": "  [!] Invalid sequence. Please try again.",
        "no_valid_nuc": "  [!] No valid nucleotide characters (A, T, C, G, U, N) found in input.",
        "step2_title": "Step 2: Select Target Organism",
        "step2_prompt": "  >> Select target organism: ",
        "invalid_organism": "  [!] Invalid organism selected. Defaulting to 'Human'.",
        "step3_title": "Step 3: Select CRISPR System",
        "step3_prompt": "  >> Select CRISPR system: ",
        "invalid_system": "  [!] Invalid system selected. Defaulting to 'SpCas9'.",
        
        # Scoring Format
        "score_ai_format": "[{model} Prediction: {score:.2f}]",
        
        "results_title": "Analysis Results",
        "target_genome": "  - Target Genome: ",
        "crispr_engine": "  - CRISPR Engine: ",
        "target_substrate": " | Target Substrate: ",
        "table_no": "No.",
        "table_strand": "Strand",
        "table_start": "Start",
        "table_end": "End",
        "table_pam": "PAM",
        "table_spacer": "Spacer Sequence (5'->3')",
        "table_score": "Score Column (Auto AI Model)",
        "total_candidates": "  Total Candidates: ",
        "scan_success": "\n  [+] Scan completed successfully. Returning to the first step...\n",
        "exiting": "\n  [!] Exiting the program. Thank you!",
        "no_candidates": "  [!] No guide RNA candidates found matching the criteria."
    },
    "ko": {
        "title_bar": "CRISPR-UNI 유전자 가위 대화형 디자인 툴",
        "welcome": "유전자 가위 디자인 플랫폼 대화형 CLI 도구에 오신 것을 환영합니다.",
        "exit_info": "어느 단계에서든 'exit'를 입력하시면 프로그램이 종료됩니다.",
        "step1_title": "1단계: 타겟 서열 입력",
        "step1_desc1": "  - 가이드를 디자인할 DNA 또는 RNA 서열을 입력해주세요.",
        "step1_desc2": "  - 예시 DNA: ATGCGATCGATCGATCGATCGATCGATCGAATCGATCGATCGATCGGGCATCGATCG",
        "step1_desc3": "  - 예시 RNA: AUCGAUCGAUCGAUCGAUCGAUCGAUGGC",
        "step1_prompt": "  >> 타겟 서열 입력: ",
        "invalid_seq": "  [!] 입력 서열이 유효하지 않습니다. 다시 시도해 주세요.",
        "no_valid_nuc": "  [!] 유효한 염기서열(A, T, C, G, U, N)이 포함되어 있지 않습니다.",
        "step2_title": "2단계: 대상 생물종 선택",
        "step2_prompt": "  >> 생물종 선택: ",
        "invalid_organism": "  [!] 올바른 생물종을 선택하지 못했습니다. 기본값인 '인간 (Human)'으로 설정합니다.",
        "step3_title": "3단계: 유전자 가위 선택",
        "step3_prompt": "  >> 유전자 가위 선택: ",
        "invalid_system": "  [!] 올바른 유전자 가위를 선택하지 못했습니다. 기본값인 'SpCas9'으로 설정합니다.",
        
        # Scoring Format
        "score_ai_format": "[{model} 모델 예측 점수: {score:.2f}]",
        
        "results_title": "분석 결과",
        "target_genome": "  - 대상 유전체: ",
        "crispr_engine": "  - 분석 도구: ",
        "target_substrate": " | 대상 substrate: ",
        "table_no": "번호",
        "table_strand": "가닥",
        "table_start": "시작",
        "table_end": "끝",
        "table_pam": "PAM/PFS",
        "table_spacer": "스페이서 서열 (5'->3')",
        "table_score": "점수 항목 (자동 AI 모델)",
        "total_candidates": "  총 후보 개수: ",
        "scan_success": "\n  [+] 스캔이 성공적으로 완료되었습니다. 첫 번째 단계로 돌아갑니다.\n",
        "exiting": "\n  [!] 프로그램을 종료합니다. 감사합니다!",
        "no_candidates": "  [!] 해당 조건에 부합하는 가이드 RNA 후보를 찾지 못했습니다."
    },
    "ja": {
        "title_bar": "CRISPR-UNI ゲノム編集対話型ツール",
        "welcome": "ゲノム編集デザインプラットフォーム対話型CLIツールへようこそ。",
        "exit_info": "どの段階でも 'exit' を入力するとプログラムが終了します。",
        "step1_title": "ステップ 1: ターゲット配列의 입력",
        "step1_desc1": "  - ガイドをデザインするDNAまたはRNAのターゲット配列を入力してください。",
        "step1_desc2": "  - DNA例: ATGCGATCGATCGATCGATCGATCGATCGAATCGATCGATCGATCGGGCATCGATCG",
        "step1_desc3": "  - RNA例: AUCGAUCGAUCGAUCGAUCGAUCGAUGGC",
        "step1_prompt": "  >> ターゲット配列入力: ",
        "invalid_seq": "  [!] 入力された配列が無効です。もう一度お試しください。",
        "no_valid_nuc": "  [!] 有効な塩基配列(A, T, C, G, U, N)が見つかりません。",
        "step2_title": "ステップ 2: 対象生物の選択",
        "step2_prompt": "  >> 対象生物を選択: ",
        "invalid_organism": "  [!] 正しい生物種が選択されませんでした。デフォルトの「ヒト (Human)」に設定します。",
        "step3_title": "ステップ 3: CRISPRシステムの選択",
        "step3_prompt": "  >> CRISPRシステムを選択: ",
        "invalid_system": "  [!] 正しいゲノム編集システムが選択されませんでした。デフォルトの「SpCas9」に設定します。",
        
        # Scoring Format
        "score_ai_format": "[{model} 予測スコア: {score:.2f}]",
        
        "results_title": "分析結果",
        "target_genome": "  - 対象ゲノム: ",
        "crispr_engine": "  - 分析エンジン: ",
        "target_substrate": " | 対象基質: ",
        "table_no": "番号",
        "table_strand": "鎖",
        "table_start": "開始",
        "table_end": "終了",
        "table_pam": "PAM/PFS",
        "table_spacer": "スペーサー配列 (5'->3')",
        "table_score": "スコア項目 (自動 AI モデル)",
        "total_candidates": "  候補総数: ",
        "scan_success": "\n  [+] スキャンが正常に完了しました。最初のステップに戻ります。\n",
        "exiting": "\n  [!] プログラムを終了します。ありがとうございました！",
        "no_candidates": "  [!] 条件に一致するガイドRNA候補が見つかりませんでした。"
    }
}

def print_separator(title: str = ""):
    if title:
        print(f"\n{'='*25} {title} {'='*25}")
    else:
        print(f"\n{'='*70}")

def print_candidates_table(candidates, lang: str = "en", predictor: BasePredictor = None):
    strings = LOCALIZATION[lang]
    if not candidates:
        print(f"  {strings['no_candidates']}")
        return
        
    no = strings["table_no"]
    strand = strings["table_strand"]
    start = strings["table_start"]
    end = strings["table_end"]
    pam = strings["table_pam"]
    spacer = strings["table_spacer"]
    score_col = strings["table_score"]
    
    print(f"  {no:<4} {strand:<6} {start:<6} {end:<5} {pam:<8} {spacer:<25} {score_col:<35}")
    print(f"  {'-'*95}")
    for idx, cand in enumerate(candidates, 1):
        if cand['score'] == "N/A":
            score_str = f"[{predictor.model_name} Score: N/A]"
        else:
            score_str = strings["score_ai_format"].format(model=predictor.model_name, score=cand['score'])
            
        print(f"  {idx:<4} {cand['strand']:<6} {cand['start']:<6} {cand['end']:<5} {cand['pam']:<8} {cand['spacer']:<25} {score_str:<35}")
    print(f"  {'-'*95}")
    print(f"{strings['total_candidates']}{len(candidates)}")

def main():
    # Configure logging here at the application entry point (not inside library modules like
    # organisms.py) so importing crispr_designer as a library never silently overrides a host
    # application's own logging configuration.
    logging.basicConfig(level=logging.INFO)

    # Print the premium Gemini-style pixel logo for CRISPR-UNI
    print(r"""
       ▄█
       ███▄       ██████╗██████╗ ██╗███████╗██████╗ ██████╗      ██╗   ██╗███╗   ██╗██║
       █████▄     ██╔════╝██╔══██╗██║██╔════╝██╔══██╗██╔══██╗     ██║   ██║████╗  ██║██║
       ███████▄   ██║     ██████╔╝██║███████╗██████╔╝██████╔╝ ══  ██║   ██║██╔██╗ ██║██║
       ███████▀   ██║     ██╔══██╗██║╚════██║██╔═══╝ ██╔══██╗     ██║   ██║██║╚██╗██║██║
       █████▀     ██████╗ ██║  ██║██║███████║██║     ██║  ██║     ╚██████╔╝██║ ╚████║██║
       ███▀       ╚═════╝ ╚═╝  ╚═╝╚═╝╚══════╝╚═╝     ╚═╝  ╚═╝      ╚═════╝ ╚═╝  ╚═══╝╚═╝
       ▀█
            [ Universal CRISPR Guide RNA Design Platform - CRISPR-UNI ]
    """)

    # Language selection prompt at startup
    print_separator("Language Selection / 언어 선택 / 言語選択")
    print("  Please select your language / 언어를 선택해 주세요 / 言語を選択してください:")
    print("  [1] English")
    print("  [2] 한국어 (Korean)")
    print("  [3] 日本語 (Japanese)")
    print()
    
    lang_choice = input("  >> Select Language (1-3): ").strip()
    if lang_choice.lower() == "exit":
        print("\n  [!] Exiting...")
        return
        
    lang = "en"
    lang_choice_lower = lang_choice.lower()
    if lang_choice == "2" or "ko" in lang_choice_lower or "kor" in lang_choice_lower or "한국" in lang_choice_lower:
        lang = "ko"
    elif lang_choice == "3" or "ja" in lang_choice_lower or "jp" in lang_choice_lower or "日本" in lang_choice_lower:
        lang = "ja"

    strings = LOCALIZATION[lang]
    print_separator(strings["title_bar"])
    print(f"  {strings['welcome']}")
    print(f"  {strings['exit_info']}\n")
    
    designer = CRISPRDesigner()
    
    # Deduplicated organism mapping for interactive selection
    organisms_list = [
        ("Human", ["Human"]),
        ("Rice", ["Rice", "벼"]),
        ("Arabidopsis", ["Arabidopsis", "애기장대"]),
        ("Mouse", ["Mouse", "마우스"]),
        ("Yeast (S. cerevisiae)", ["Yeast", "Saccharomyces cerevisiae", "효모"]),
        ("E. coli", ["E. coli", "Ecoli", "대장균"])
    ]
    
    # Configured CRISPR systems in registry
    systems_list = list(CRISPR_REGISTRY.keys())
    
    while True:
        # Step 1: Enter Target Sequence
        print_separator(strings["step1_title"])
        print(strings["step1_desc1"])
        print(strings["step1_desc2"])
        print(strings["step1_desc3"] + "\n")
        
        sequence_input = input(strings["step1_prompt"]).strip()
        if sequence_input.lower() == "exit":
            print(strings["exiting"])
            break
            
        if not sequence_input:
            print(strings["invalid_seq"])
            continue
            
        # Clean sequence characters
        cleaned_sequence = "".join(char for char in sequence_input if char.upper() in 'ATCGUN')
        if not cleaned_sequence:
            print(strings["no_valid_nuc"])
            continue
            
        # Step 2: Select Target Organism
        print_separator(strings["step2_title"])
        for idx, (display_name, _) in enumerate(organisms_list, 1):
            print(f"  [{idx}] {display_name}")
        print()
        
        organism_choice = input(strings["step2_prompt"]).strip()
        if organism_choice.lower() == "exit":
            print(strings["exiting"])
            break
            
        # Resolve organism selection
        selected_organism = None
        # Try matching index
        if organism_choice.isdigit():
            idx = int(organism_choice) - 1
            if 0 <= idx < len(organisms_list):
                selected_organism = organisms_list[idx][1][0]
        else:
            # Try matching name strings
            for display, aliases in organisms_list:
                if any(alias.lower() in organism_choice.lower() for alias in aliases):
                    selected_organism = aliases[0]
                    break
        
        if not selected_organism:
            print(strings["invalid_organism"])
            selected_organism = "Human"
            
        # Step 3: Select CRISPR System
        print_separator(strings["step3_title"])
        for idx, sys_name in enumerate(systems_list, 1):
            sys_instance = CRISPR_REGISTRY[sys_name]
            target_desc = getattr(sys_instance, 'target_type', 'DNA')
            print(f"  [{idx}] {sys_name:<15} (PAM/PFS: {sys_instance.pam:<8} | Spacer: {sys_instance.spacer_length}nt | Target: {target_desc})")
        print()
        
        system_choice = input(strings["step3_prompt"]).strip()
        if system_choice.lower() == "exit":
            print(strings["exiting"])
            break
            
        # Resolve system selection
        selected_system = None
        if system_choice.isdigit():
            idx = int(system_choice) - 1
            if 0 <= idx < len(systems_list):
                selected_system = systems_list[idx]
        else:
            for sys_name in systems_list:
                if sys_name.lower() in system_choice.lower():
                    selected_system = sys_name
                    break
                    
        if not selected_system:
            print(strings["invalid_system"])
            selected_system = "SpCas9"

        # Apply selection to designer and search
        designer.set_organism(selected_organism)
        designer.set_crispr_system(selected_system)
        
        # Get active system instance and its dynamic deep learning predictor
        sys_instance = CRISPR_REGISTRY[selected_system]
        predictor = PredictorFactory.get_predictor(selected_system)
        
        # Perform Guide RNA Scan
        candidates = designer.find_candidates(cleaned_sequence)
        
        # Polymorphically predict efficiency for each guide RNA candidate
        weights_missing_notified = False
        for cand in candidates:
            kwargs = {}
            if selected_system == 'Prime_Editor':
                # Prime editor needs PBS and RT template length bounds
                kwargs['pbs_len'] = getattr(sys_instance, 'pbs_length', 13)
                kwargs['rt_len'] = getattr(sys_instance, 'rt_length', 15)
            elif selected_system == 'Cas13d':
                # Cas13d RNA models evaluates Protospacer Flanking Site (PFS)
                kwargs['PFS'] = cand['pam']
            else:
                # DNA models evaluates standard seed double-strand break (DSB) parameters
                kwargs['target_site'] = cand['target_site']
                
            try:
                # Calculate scores via the concrete deep learning model
                cand["score"] = predictor.predict_efficiency(spacer=cand['spacer'], **kwargs)
            except Exception as e:
                if isinstance(e, ModelWeightsMissingError):
                    if not weights_missing_notified:
                        file_name = getattr(e, 'missing_file', predictor.weight_file)
                        dl_url = getattr(e, 'download_url', predictor.download_url)
                        if lang == "ko":
                            print(f"\n  ❌ [{predictor.model_name}] 가중치 파일이 없습니다. 경로를 확인하세요.")
                            print(f"      다운로드 링크: {dl_url}\n")
                        elif lang == "ja":
                            print(f"\n  ❌ [{predictor.model_name}] 重みファイルがありません。パスを確認してください。")
                            print(f"      ダウンロードリンク: {dl_url}\n")
                        else:
                            print(f"\n  ❌ [{predictor.model_name}] Weights file is missing. Please verify path.")
                            print(f"      Download Link: {dl_url}\n")
                        weights_missing_notified = True
                cand["score"] = "N/A"
        
        # Display analysis progress
        print_separator(strings["results_title"])
        print(f"{strings['target_genome']}{designer.organism_name} ({designer.organism_metadata['ref_assembly']})")
        print(f"{strings['crispr_engine']}{designer.active_system.name} (PAM/PFS: {designer.active_system.pam}{strings['target_substrate']}{getattr(designer.active_system, 'target_type', 'DNA')})")
        print()
        
        # Print table
        print_candidates_table(candidates, lang, predictor)
        print(strings["scan_success"])

if __name__ == "__main__":
    main()
