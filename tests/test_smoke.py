import aniso_brem


def test_eebls_returns_positive_float():
    value = aniso_brem.eebls(1500.0, 100.0, 200)
    assert isinstance(value, float)
    assert value > 0.0


def test_cspe_returns_positive_float():
    value = aniso_brem.cspe(1500.0, 100.0)
    assert isinstance(value, float)
    assert value > 0.0
