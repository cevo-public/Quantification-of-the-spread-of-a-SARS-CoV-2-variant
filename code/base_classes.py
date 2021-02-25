from dataclasses import dataclass
from typing import List, Any

import pandas as pd


@dataclass
class Dataset:
    # columns: date:datetime, year:int, week:int, t:int, region:str (optional), cases:int
    cases: pd.DataFrame

    # columns: date:datetime, region:str (optional), sequenced:int, b117:int, original:int
    variants: pd.DataFrame


@dataclass
class PlotConfig:
    colors: List[str]
    greys: List[str]


@dataclass(frozen=True)
class ConfidenceInterval:
    lower: float
    upper: float


@dataclass(frozen=True)
class GrowthNumbers:
    alpha: float
    generation_time: float
    reproduction_number: float
    a_mle: float
    a_ci: ConfidenceInterval
    t0_mle: float
    t0_ci: ConfidenceInterval
    fd_mle: float  # Discrete setting, exp(a*g)-1
    fd_ci: ConfidenceInterval
    fc_mle: float  # Continuous setting, a*g
    fc_ci: ConfidenceInterval
    statsmodel_model: Any
