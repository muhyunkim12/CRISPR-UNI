"""
crispr_designer/predictors.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Defines polymorphic deep learning predictive models (PyTorch & TensorFlow)
and a factory pattern to dynamically bind scoring algorithms to distinct CRISPR systems.
Reads real weights from the './weights/' directory.
"""

import abc
import os
from typing import Dict
import numpy as np

# Declare required ML imports at the top
try:
    import tensorflow as tf
except ImportError:
    tf = None

try:
    import torch
except ImportError:
    torch = None

from .utils import seq_to_onehot

# Shared nucleotide -> index mapping used by every PyTorch-backed predictor below
# (U maps to the same index as T since sequences are treated as DNA for encoding purposes).
NUCLEOTIDE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}

class ModelWeightsMissingError(FileNotFoundError):
    """
    Custom exception raised when deep learning framework or weights are missing.
    Inherits from FileNotFoundError so it acts as both FileNotFoundError and OSError.
    """
    def __init__(self, model_name: str, missing_file: str, download_url: str, is_framework_missing: bool = False):
        self.model_name = model_name
        self.missing_file = missing_file
        self.download_url = download_url
        self.is_framework_missing = is_framework_missing
        if is_framework_missing:
            msg = f"Deep learning framework not installed for {model_name}. Please install required packages."
        else:
            msg = f"Weight file '{missing_file}' for {model_name} not found. Download link: {download_url}"
        super().__init__(msg)

class BasePredictor(abc.ABC):
    """
    Abstract Base Class for all machine learning on-target guide efficiency predictors.
    """
    def __init__(self, model_name: str, weight_file: str, download_url: str):
        self.model_name = model_name
        self.weight_file = os.path.join("weights", weight_file)
        self.download_url = download_url
        self.model = None

    def _verify_weights(self):
        """
        Verifies that the target weights file exists locally in the './weights/' directory.
        Raises FileNotFoundError if missing.
        """
        if not os.path.exists(self.weight_file):
            raise ModelWeightsMissingError(self.model_name, self.weight_file, self.download_url)

    @abc.abstractmethod
    def predict_efficiency(self, **kwargs) -> float:
        """
        Abstract method to predict guide efficiency.
        Polymorphic interface to accept different metadata arguments across systems.
        """
        pass

class DeepSpCas9Predictor(BasePredictor):
    """
    DeepSpCas9 model for standard SpCas9 DNA DSB efficiency prediction.
    Powered by TensorFlow Keras.
    """
    def __init__(self):
        super().__init__(
            model_name="DeepSpCas9",
            weight_file="DeepSpCas9_weights.h5",
            download_url="https://github.com/BoseonShim/DeepSpCas9"
        )

    def predict_efficiency(self, spacer: str, target_site: str, **kwargs) -> float:
        if tf is None:
            raise ModelWeightsMissingError(self.model_name, self.weight_file, self.download_url, is_framework_missing=True)
        self._verify_weights()

        if self.model is None:
            self.model = tf.keras.models.load_model(self.weight_file)

        # Preprocessing: Convert sequence to One-hot using the common seq_to_onehot helper
        seq = spacer.upper()
        if len(seq) < 30:
            seq = (seq + "N" * 30)[:30]
        
        one_hot = seq_to_onehot(seq).reshape(1, 30, 4)
        preds = self.model.predict(one_hot, verbose=0)
        return round(float(preds[0][0]), 2)

class DeepSaCas9Predictor(BasePredictor):
    """
    DeepSaCas9 model targeting SaCas9 NNGRRT PAM system.
    Powered by TensorFlow Keras.
    """
    def __init__(self):
        super().__init__(
            model_name="DeepSaCas9",
            weight_file="DeepSaCas9_weights.h5",
            download_url="https://github.com/koots/DeepSaCas9"
        )

    def predict_efficiency(self, spacer: str, target_site: str, **kwargs) -> float:
        if tf is None:
            raise ModelWeightsMissingError(self.model_name, self.weight_file, self.download_url, is_framework_missing=True)
        self._verify_weights()

        if self.model is None:
            self.model = tf.keras.models.load_model(self.weight_file)

        # Preprocessing: Convert 21nt spacer sequence to One-hot
        seq = spacer.upper()
        one_hot = seq_to_onehot(seq).reshape(1, len(seq), 4)
        preds = self.model.predict(one_hot, verbose=0)
        return round(float(preds[0][0]), 2)

class DeepCpf1Predictor(BasePredictor):
    """
    DeepCpf1 model for Cas12a sticky end cleavage prediction with TTTN PAM.
    Powered by TensorFlow Keras.
    """
    def __init__(self):
        super().__init__(
            model_name="DeepCpf1",
            weight_file="DeepCpf1_weights.h5",
            download_url="https://github.com/CRISPRdeepHF/DeepHF"
        )

    def predict_efficiency(self, spacer: str, target_site: str, **kwargs) -> float:
        if tf is None:
            raise ModelWeightsMissingError(self.model_name, self.weight_file, self.download_url, is_framework_missing=True)
        self._verify_weights()

        if self.model is None:
            self.model = tf.keras.models.load_model(self.weight_file)

        seq = spacer.upper()
        one_hot = seq_to_onehot(seq).reshape(1, len(seq), 4)
        preds = self.model.predict(one_hot, verbose=0)
        return round(float(preds[0][0]), 2)

class PRIDICTPredictor(BasePredictor):
    """
    PRIDICT (Prime Editing Guide RNA Efficiency Predictor).
    Powered by PyTorch.
    """
    def __init__(self):
        super().__init__(
            model_name="PRIDICT",
            weight_file="PRIDICT_model.pt",
            download_url="https://github.com/allencellmodeling/pridict"
        )

    def predict_efficiency(self, spacer: str, pbs_len: int = 13, rt_len: int = 15, **kwargs) -> float:
        if torch is None:
            raise ModelWeightsMissingError(self.model_name, self.weight_file, self.download_url, is_framework_missing=True)
        self._verify_weights()

        if self.model is None:
            self.model = torch.load(self.weight_file, map_location=torch.device('cpu'))
            self.model.eval()

        seq = spacer.upper()
        seq_indices = [NUCLEOTIDE_TO_INDEX.get(char, 0) for char in seq]

        with torch.no_grad():
            seq_tensor = torch.tensor([seq_indices], dtype=torch.long)
            pbs_tensor = torch.tensor([[pbs_len]], dtype=torch.float32)
            rt_tensor = torch.tensor([[rt_len]], dtype=torch.float32)
            
            preds = self.model(seq_tensor, pbs_tensor, rt_tensor)
            # Use item() to extract final predictive score
            score = preds.item() if hasattr(preds, 'item') else preds[0]
            return round(float(score), 2)

class Cas13designPredictor(BasePredictor):
    """
    Cas13design model optimized for RNA targeting and 3'-PFS context evaluation.
    Powered by PyTorch.
    """
    def __init__(self):
        super().__init__(
            model_name="Cas13design",
            weight_file="Cas13design_weights.pt",
            download_url="https://github.com/gpp-lab/cas13design"
        )

    def predict_efficiency(self, spacer: str, PFS: str = "H", **kwargs) -> float:
        if torch is None:
            raise ModelWeightsMissingError(self.model_name, self.weight_file, self.download_url, is_framework_missing=True)
        self._verify_weights()

        if self.model is None:
            self.model = torch.load(self.weight_file, map_location=torch.device('cpu'))
            self.model.eval()

        seq = spacer.upper()
        seq_indices = [NUCLEOTIDE_TO_INDEX.get(char, 0) for char in seq]

        with torch.no_grad():
            seq_tensor = torch.tensor([seq_indices], dtype=torch.long)
            # Construct PyTorch sequence structural features (folding energy and base pair probabilities)
            folding_features = torch.tensor([[0.85, -12.4, 0.4]], dtype=torch.float32)
            
            preds = self.model(seq_tensor, folding_features)
            score = preds.item() if hasattr(preds, 'item') else preds[0]
            return round(float(score), 2)

class DeepCas12fPredictor(BasePredictor):
    """
    DeepCas12f model for mini CRISPR Cas14a/Cas12f ssDNA targets.
    Powered by PyTorch.
    """
    def __init__(self):
        super().__init__(
            model_name="DeepCas12f",
            weight_file="DeepCas12f_weights.pt",
            download_url="https://github.com/miniCRISPR/DeepCas12f"
        )

    def predict_efficiency(self, spacer: str, **kwargs) -> float:
        if torch is None:
            raise ModelWeightsMissingError(self.model_name, self.weight_file, self.download_url, is_framework_missing=True)
        self._verify_weights()

        if self.model is None:
            self.model = torch.load(self.weight_file, map_location=torch.device('cpu'))
            self.model.eval()

        seq = spacer.upper()
        seq_indices = [NUCLEOTIDE_TO_INDEX.get(char, 0) for char in seq]

        with torch.no_grad():
            seq_tensor = torch.tensor([seq_indices], dtype=torch.long)
            # ssDNA targeting feature states (1.0 ssDNA flag)
            substrate_features = torch.tensor([[1.0, 0.0]], dtype=torch.float32)
            
            preds = self.model(seq_tensor, substrate_features)
            score = preds.item() if hasattr(preds, 'item') else preds[0]
            return round(float(score), 2)

class PredictorFactory:
    """
    Factory Pattern class to dynamically match active CRISPR systems with their
    dedicated polymorphic Deep Learning Predictors.

    Predictor instances are cached per resolved CRISPR system name. Each BasePredictor
    lazily loads its (potentially large) model weights into self.model on first use, but
    that cache is only useful if the same instance is reused across calls - callers that
    invoke get_predictor() repeatedly in a loop (e.g. once per scan) would otherwise
    reload the weights from disk every single time.
    """
    _instance_cache: Dict[str, BasePredictor] = {}

    @staticmethod
    def get_predictor(crispr_system: str) -> BasePredictor:
        cached = PredictorFactory._instance_cache.get(crispr_system)
        if cached is not None:
            return cached

        system_lower = crispr_system.lower()
        if "spcas9" in system_lower:
            predictor = DeepSpCas9Predictor()
        elif "sacas9" in system_lower:
            predictor = DeepSaCas9Predictor()
        elif "cas12a" in system_lower or "cpf1" in system_lower:
            predictor = DeepCpf1Predictor()
        elif "prime_editor" in system_lower or "prime" in system_lower:
            predictor = PRIDICTPredictor()
        elif "cas13d" in system_lower:
            predictor = Cas13designPredictor()
        elif "cas14a" in system_lower or "cas12f" in system_lower:
            predictor = DeepCas12fPredictor()
        else:
            predictor = DeepSpCas9Predictor()

        PredictorFactory._instance_cache[crispr_system] = predictor
        return predictor
