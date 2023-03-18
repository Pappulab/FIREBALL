"""
Unit and regression test for the fireball package.
"""

# Import package, test suite, and other packages as needed
import fireball
import pytest
import sys

def test_fireball_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "fireball" in sys.modules
