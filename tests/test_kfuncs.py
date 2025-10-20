import pytest
from src.bb_ions.kfuncs import *


@pytest.mark.parametrize("l, m", [
    (3, 5),
    (12, 6),
    (12, 12),
])

def test_convtok(l, m):
    k = 0
    for u in range(4):
        for v in range(l):
            for w in range(m):
                assert(convtok(l, m, u, v, w) == k)
                k = k + 1