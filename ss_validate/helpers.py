import math
import pandas as pd
import numpy as np
from pandas_schema.validation import (MatchesPatternValidation,
                                      InListValidation,
                                      CanConvertValidation,
                                      CustomSeriesValidation,
                                      _SeriesValidation)
from ss_validate import __version__


class InInclusiveRangeValidation(_SeriesValidation):
    """
    Checks that each element in the series is within a given inclusive
    numerical range like so x <= n < y.
    Doesn't care if the values are not numeric - it will try anyway.
    """
    def __init__(self, min: float = -math.inf, max: float = math.inf, **kwargs):
        """
        :param min: The minimum (inclusive) value to accept
        :param max: The maximum (inclusive) value to accept
        """
        self.min = min
        self.max = max
        super().__init__(**kwargs)

    @property
    def default_message(self):
        return 'was not >= {} and <= {})'.format(self.min, self.max)

    def validate(self, series: pd.Series) -> pd.Series:
        series = pd.to_numeric(series, errors='coerce')
        return (series >= self.min) & (series <= self.max)


class InExclusiveRangeValidation(_SeriesValidation):
    """
    Checks that each element in the series is within a given exclusive numerical range.
    """
    def __init__(self, min: float = -math.inf, max: float = math.inf, **kwargs):
        """
        :param min: The minimum (exclusive) value to accept
        :param max: The maximum (exclusive) value to accept
        """
        self.min = min
        self.max = max
        super().__init__(**kwargs)

    @property
    def default_message(self):
        return 'was not > {} and < {})'.format(self.min, self.max)

    def validate(self, series: pd.Series) -> pd.Series:
        series = pd.to_numeric(series, errors='coerce')
        return (series > self.min) & (series < self.max)


class InRangeValidationUpperInclusive(_SeriesValidation):
    """
    Checks that each element in the series is within a given numerical range.
    """
    def __init__(self, min: float = -math.inf, max: float = math.inf, **kwargs):
        """
        :param min: The minimum (exclusive) value to accept
        :param max: The maximum (inclusive) value to accept
        """
        self.min = min
        self.max = max
        super().__init__(**kwargs)

    @property
    def default_message(self):
        return 'was not > {} and <= {})'.format(self.min, self.max)

    def validate(self, series: pd.Series) -> pd.Series:
        series = pd.to_numeric(series, errors='coerce')
        return (series > self.min) & (series <= self.max)


def match_regex(pattern):
    return MatchesPatternValidation(r'{}'.format(pattern))


def in_list(list):
    return InListValidation(list, message=f"is not an accepted value. Value must be be {list}")


def in_range(lower=-math.inf, upper=math.inf):
    return InInclusiveRangeValidation(lower, upper)


def is_dtype(dtype):
    return CanConvertValidation(dtype)


p_value_validation = InRangeValidationUpperInclusive(0, 1) | (
        CustomSeriesValidation(
            lambda x: pd.to_numeric(x.str.split('e|E', expand=True)[1].fillna(value=np.nan)
                                    , errors='coerce') < -1,
            'Numbers should be between 0 and 1') &
        CustomSeriesValidation(
            lambda x: pd.to_numeric(x.str.split('e|E', expand=True)[0].fillna(value=np.nan)
                                    , errors='coerce') > 0,
            'Numbers should be between 0 and 1')
)

p_value_validation_allow_zero = InInclusiveRangeValidation(0, 1) | (
        CustomSeriesValidation(
            lambda x: pd.to_numeric(x.str.split('e|E', expand=True)[1].fillna(value=np.nan)
                                    , errors='coerce') < -1,
            'Numbers should be between 0 and 1') &
        CustomSeriesValidation(
            lambda x: pd.to_numeric(x.str.split('e|E', expand=True)[0].fillna(value=np.nan)
                                    , errors='coerce') > 0,
            'Numbers should be between 0 and 1')
)


def get_version():
    return __version__
