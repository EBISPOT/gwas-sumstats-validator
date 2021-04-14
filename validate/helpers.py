import math
import pandas as pd
from pandas_schema.validation import _SeriesValidation


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
