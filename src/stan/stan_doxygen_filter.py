#!/usr/bin/env python3

import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")

# Convert Stan tuple return types to C++-style template types.
#
# Example:
#   tuple(
#     vector,
#     real,
#     matrix
#   )
#
# becomes:
#   std::tuple<
#     vector,
#     real,
#     matrix
#   >
#
# The pattern is restricted to tuple contents without parentheses, which
# matches Stan type lists such as vector, real, matrix, and int.
text = re.sub(
    r"\btuple\s*\(([^()]*)\)",
    r"std::tuple<\1>",
    text,
    flags=re.DOTALL,
)

sys.stdout.write(text)
