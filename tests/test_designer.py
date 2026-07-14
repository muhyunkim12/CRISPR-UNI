"""
tests/test_designer.py
~~~~~~~~~~~~~~~~~~~~~~
Automated unit test suite for the Universal CRISPR Guide RNA Design Platform.
Validates bioinformatic helpers, dynamic registry behavior, coordinate mapping, and runtime extensibility.
"""

import unittest
import os
from crispr_designer import (
    CRISPRDesigner,
    CRISPRSystem,
    CRISPR_REGISTRY,
    ORGANISM_REGISTRY
)
from crispr_designer.utils import reverse_complement, iupac_to_regex
from crispr_designer.systems import SpCas9System, Cas12aSystem, SaCas9System

class TestUtils(unittest.TestCase):
    """
    Tests bioinformatic helper utilities.
    """
    def test_reverse_complement(self):
        # Standard DNA base test
        self.assertEqual(reverse_complement("ATCG"), "CGAT")
        # Case sensitivity test
        self.assertEqual(reverse_complement("aTcG"), "CgAt")
        # Degenerate base test (R maps to Y, etc.)
        self.assertEqual(reverse_complement("ATGRY"), "RYCAT")
        # Dual strand identity test
        self.assertEqual(reverse_complement(reverse_complement("ATCG")), "ATCG")

    def test_iupac_to_regex(self):
        # Precise bases
        self.assertEqual(iupac_to_regex("GG"), "GG")
        # Degenerates
        self.assertEqual(iupac_to_regex("NGG"), "[ACGT]GG")
        # Multi-degenerate
        self.assertEqual(iupac_to_regex("NNGRRT"), "[ACGT][ACGT]G[AG][AG]T")


class TestCRISPRSystems(unittest.TestCase):
    """
    Tests for individual CRISPR Endonuclease system engines.
    """
    def test_spcas9_properties(self):
        system = SpCas9System()
        self.assertEqual(system.name, "SpCas9")
        self.assertEqual(system.spacer_length, 20)
        self.assertEqual(system.pam_position, "3_prime")

    def test_spcas9_scanning(self):
        system = SpCas9System()
        # forward target sequence: 20nt spacer + PAM (NGG)
        seq = "ATCGATCGATCGATCGATCGGGG"  # 'ATCGATCGATCGATCGATCG' is 20nt, 'GGG' matches NGG
        candidates = system.find_candidates(seq)
        
        self.assertGreaterEqual(len(candidates), 1)
        # Check candidate metadata mapping
        cand = candidates[0]
        self.assertEqual(cand['strand'], '+')
        self.assertEqual(cand['spacer'], "ATCGATCGATCGATCGATCG")
        self.assertEqual(cand['pam'], "GGG")
        self.assertEqual(cand['start'], 0)
        self.assertEqual(cand['end'], 23)

    def test_azimuth_context_extraction(self):
        system = SpCas9System()
        # Upstream 4nt: AAAA, Spacer 20nt: ATCGATCGATCGATCGATCG, PAM 3nt: TGG, Downstream 3nt: TTT
        seq = "AAAAATCGATCGATCGATCGATCGTGGTTT"
        candidates = system.find_candidates(seq)
        self.assertEqual(len(candidates), 1)
        cand = candidates[0]
        self.assertEqual(cand['context_seq_30'], "AAAAATCGATCGATCGATCGATCGTGGTTT")

        # Test boundary scenario where sequence lacks upstream 4nt (only 3nt)
        seq_short = "AAATCGATCGATCGATCGATCGTGGTTT"
        candidates_short = system.find_candidates(seq_short)
        self.assertEqual(len(candidates_short), 1)
        self.assertIsNone(candidates_short[0]['context_seq_30'])

    def test_cas12a_scanning(self):
        system = Cas12aSystem()
        self.assertEqual(system.pam_position, "5_prime")
        
        # 5' PAM: TTTN + 23nt Spacer
        # seq has 'TTTA' (4nt) + 'ATCGATCGATCGATCGATCGATC' (23nt)
        seq = "TTTAATCGATCGATCGATCGATCGATC"
        candidates = system.find_candidates(seq)
        self.assertEqual(len(candidates), 1)
        
        cand = candidates[0]
        self.assertEqual(cand['strand'], '+')
        self.assertEqual(cand['pam'], "TTTA")
        self.assertEqual(cand['spacer'], "ATCGATCGATCGATCGATCGATC")
        self.assertEqual(cand['start'], 0)
        self.assertEqual(cand['end'], 27)

    def test_cas14a_scanning(self):
        from crispr_designer import Cas14aSystem
        system = Cas14aSystem()
        self.assertEqual(system.target_type, "ssDNA")
        
        # 5' PAM: TTTA + 20nt Spacer
        seq = "TTTAATCGATCGATCGATCGATCG"
        candidates = system.find_candidates(seq)
        self.assertEqual(len(candidates), 1)
        cand = candidates[0]
        self.assertEqual(cand['spacer'], "ATCGATCGATCGATCGATCG")
        self.assertEqual(cand['pam'], "TTTA")
        self.assertEqual(cand['start'], 0)
        self.assertEqual(cand['end'], 24)


class TestOrganisms(unittest.TestCase):
    """
    Tests reference genome settings and NCBI accession conversions.
    """
    def test_organism_mapping(self):
        # Verify Korean mapping key resolves to English equivalents
        self.assertIn("인간", ORGANISM_REGISTRY)
        self.assertIn("Human", ORGANISM_REGISTRY)
        
        human_meta = ORGANISM_REGISTRY["인간"]
        self.assertEqual(human_meta["ref_assembly"], "GRCh38")
        
        rice_meta = ORGANISM_REGISTRY["벼"]
        self.assertEqual(rice_meta["ref_assembly"], "IRGSP-1.0")

    def test_new_organisms(self):
        # Verify GRCm39, sacCer3, K-12 mappings
        self.assertEqual(ORGANISM_REGISTRY["마우스"]["ref_assembly"], "GRCm39")
        self.assertEqual(ORGANISM_REGISTRY["Mouse"]["ref_assembly"], "GRCm39")
        self.assertEqual(ORGANISM_REGISTRY["효모"]["ref_assembly"], "sacCer3")
        self.assertEqual(ORGANISM_REGISTRY["Saccharomyces cerevisiae"]["ref_assembly"], "sacCer3")
        self.assertEqual(ORGANISM_REGISTRY["대장균"]["ref_assembly"], "K-12")
        self.assertEqual(ORGANISM_REGISTRY["E. coli"]["ref_assembly"], "K-12")


class TestCRISPRDesigner(unittest.TestCase):
    """
    Tests main orchestrator CRISPRDesigner and runtime extensibility.
    """
    def test_designer_lifecycle(self):
        designer = CRISPRDesigner()
        
        # Ensure setting organisms and systems works
        designer.set_organism("벼")
        self.assertEqual(designer.organism_name, "벼")
        
        designer.set_crispr_system("Cas12a")
        self.assertEqual(designer.active_system.name, "Cas12a")
        
    def test_custom_runtime_system_registration(self):
        designer = CRISPRDesigner()
        
        # Define a custom Cas-X at runtime
        class CustomCasX(CRISPRSystem):
            def __init__(self):
                super().__init__(
                    name="CustomCasX",
                    pam="CCA",
                    spacer_length=15,
                    pam_position="3_prime"
                )
        
        # Register the custom system dynamically
        custom_instance = CustomCasX()
        designer.register_custom_system("CustomCasX", custom_instance)
        
        # Verify it is registered in global database
        self.assertIn("CustomCasX", CRISPR_REGISTRY)
        self.assertEqual(designer.active_system.name, "CustomCasX")
        
        # Scan with it
        test_seq = "ATCGATCGATCGATCCCAGGG" # 15nt spacer + 'CCA' PAM
        candidates = designer.find_candidates(test_seq)
        self.assertGreaterEqual(len(candidates), 1)
        
        cand = candidates[0]
        self.assertEqual(cand['pam'], "CCA")
        self.assertEqual(cand['spacer'], "ATCGATCGATCGATC")


class TestDeepPredictors(unittest.TestCase):
    """
    Tests the new PredictorFactory and deep learning predictors architecture.
    """
    def test_factory_matching(self):
        from crispr_designer import PredictorFactory
        
        # Test SpCas9 matching
        pred_sp = PredictorFactory.get_predictor("SpCas9")
        self.assertEqual(pred_sp.model_name, "DeepSpCas9")
        
        # Test Cas12a matching
        pred_cpf1 = PredictorFactory.get_predictor("Cas12a")
        self.assertEqual(pred_cpf1.model_name, "DeepCpf1")

        # Test Prime Editor matching
        pred_pri = PredictorFactory.get_predictor("Prime_Editor")
        self.assertEqual(pred_pri.model_name, "PRIDICT")

        # SaCas9 and Cas14a(Cas12f) have no real, publicly available pretrained model -
        # the factory returns None rather than guessing with an unrelated model.
        self.assertIsNone(PredictorFactory.get_predictor("SaCas9"))
        self.assertIsNone(PredictorFactory.get_predictor("Cas14a(Cas12f)"))

    def test_polymorphic_prediction_scores(self):
        from crispr_designer import PredictorFactory
        from crispr_designer.predictors import ModelWeightsMissingError
        
        # Test SpCas9 predict call raises ModelWeightsMissingError due to missing weights
        pred_sp = PredictorFactory.get_predictor("SpCas9")
        with self.assertRaises(ModelWeightsMissingError) as ctx:
            pred_sp.predict_efficiency(spacer="ATCGATCGATCGATCGATCG", target_site="ATCGATCGATCGATCGATCGTGG")
        self.assertEqual(ctx.exception.missing_file, os.path.join("weights", "DeepSpCas9_weights.h5"))
        self.assertEqual(ctx.exception.download_url, "https://github.com/BoseonShim/DeepSpCas9")
        
        # Test PRIDICT Prime Editor raises ModelWeightsMissingError due to missing weights
        pred_pri = PredictorFactory.get_predictor("Prime_Editor")
        with self.assertRaises(ModelWeightsMissingError) as ctx:
            pred_pri.predict_efficiency(spacer="ATCGATCGATCGATCGATCG", pbs_len=13, rt_len=15)
        self.assertEqual(ctx.exception.missing_file, os.path.join("weights", "PRIDICT_model.pt"))
        self.assertEqual(ctx.exception.download_url, "https://github.com/allencellmodeling/pridict")


class TestSaCas9SpecificityScore(unittest.TestCase):
    """
    Tests calc_sacas9_score against known values from the real, published SaCas9
    specificity model table (Tycko et al. 2018; see crispr_designer/sacas9_scores.py).
    """
    def test_perfect_match_scores_one(self):
        from crispr_designer import calc_sacas9_score
        guide = "A" * 21
        self.assertEqual(calc_sacas9_score(guide, guide), 1.0)

    def test_position_one_mismatch_is_ignored(self):
        from crispr_designer import calc_sacas9_score
        guide = "A" * 21
        # Mismatch only at the 5'-most position (position 1) - the model explicitly does not
        # score this position, so it should be treated as if there were no mismatch at all.
        hit = "C" + "A" * 20
        self.assertEqual(calc_sacas9_score(guide, hit), 1.0)

    def test_single_mismatch_matches_table_value(self):
        from crispr_designer import calc_sacas9_score
        from crispr_designer.sacas9_scores import SACAS9_MISMATCH_SCORES
        guide = "A" * 21
        # Mismatch at position 5 (1-indexed): hit has 'C' where guide has 'A'.
        hit = "A" * 4 + "C" + "A" * 16
        expected = SACAS9_MISMATCH_SCORES['CA,5']
        self.assertAlmostEqual(calc_sacas9_score(guide, hit), expected)

    def test_wrong_length_returns_none(self):
        from crispr_designer import calc_sacas9_score
        self.assertIsNone(calc_sacas9_score("A" * 20, "A" * 20))
        self.assertIsNone(calc_sacas9_score("A" * 21, "A" * 20))


if __name__ == "__main__":
    unittest.main()
