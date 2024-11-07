"""Top-level package for gdb_brexazine."""

__author__ = """Ye Buehler"""
__email__ = 'ye.buehler-feng@unibe.ch'
__version__ = '0.1.0'


# Add imports here
from .utils import Utils

# When using <<from gdb_brexazine import * >> only the classes listed in __all__ will be imported.
__all__ = [
    'Utils',
]

