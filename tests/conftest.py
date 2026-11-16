"""Hypothesis profiles for the fuzz tests.

* default: quick enough for every test run
* fuzz: many more examples, for the dedicated CI job or a local soak test:
  ``pytest tests/test_fuzz.py --hypothesis-profile=fuzz``
"""

from hypothesis import HealthCheck, settings

# deadline=None: an example's running time varies with sequence length, which is not a bug
settings.register_profile("default", max_examples=300, deadline=None)
settings.register_profile("fuzz", max_examples=20_000, deadline=None, suppress_health_check=[HealthCheck.too_slow])
settings.load_profile("default")
