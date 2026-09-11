"""its-fun-2-map: a pipeline for extracting and validating fungal ITS barcode
sequences from museum specimen genome skims.

Each pipeline step lives in its own module and is exposed as an ``itsfun-*``
console script. Nothing is re-exported here: modules are imported by their own
path (e.g. ``from its_fun_2_map import fastp_module``) so that importing the
package does not pull in heavy optional dependencies such as rpy2.
"""

import logging

__version__ = "1.0.0"

# Library-safe logging: adding a NullHandler means importing this package has no
# logging side effects. Each module's main() calls setup_logging() to configure
# handlers for CLI use only.
logging.getLogger(__name__).addHandler(logging.NullHandler())
