#!/usr/bin/env python3

import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")

# Present Stan tuple types as a simple return type for Doxygen's C++ parser.
#
# Example:
#   tuple(
#     vector,
#     real,
#     matrix
#   )
#
# becomes:
#   tuple
#
# The pattern is restricted to tuple contents without parentheses, which
# matches Stan type lists such as vector, real, matrix, and int.
# Retain line numbers so source links still point to the original Stan lines.
text = re.sub(
    r"\btuple\s*\(([^()]*)\)",
    lambda match: "tuple" + "\n" * match.group(0).count("\n"),
    text,
    flags=re.DOTALL,
)

sys.stdout.write(text)
