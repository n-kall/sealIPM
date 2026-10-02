# Building the Stan HTML reference

From the repository root:

```sh
cd src/stan
doxygen Doxyfile
```

Open `html/index.html`. Rebuild after changing function documentation.
Doxygen 1.9.8 or newer, Python 3, and Graphviz (`dot`) are required.
The paths in `Doxyfile` are relative to `src/stan`, so run the command there.

The configuration follows the Stan-as-C++ parsing, search, source browsing,
tree navigation, and MathJax approach used in
[epinowcast's Doxyfile](https://github.com/epinowcast/epinowcast/blob/main/inst/stan/Doxyfile).
The separate overview page, layout, and stylesheet also draw on the organization
described in the epiforecasts guide
[How we document Stan functions](https://epiforecasts.io/posts/2025-03-19-stan-doc-guide/).
The local stylesheet is original; no external theme is needed at build time.

MathJax is pinned to version 3.2.2 on a CDN. To display equations offline,
install MathJax locally and update `MATHJAX_RELPATH` to its package directory
relative to the generated HTML directory.

Edit `mainpage.md`, `DoxygenLayout.xml`, or `sealipm.css` here, not the generated
files in `html/`. The documentation filter affects parsing only: source
listings retain the original Stan syntax. No model equations are changed.

Tuple return types are displayed as `tuple` rather than a C++ template signature
The `@return` description gives the contents and their order
The filter preserves newlines so source links retain the original line numbers

Doxygen joins consecutive prose lines into a single paragraph
Use an empty comment line (` *`) between separate statements or paragraphs
Keep wrapped sentences together and use `<br>` only when a line break within
the same paragraph is intentional
See [Doxygen paragraph formatting](https://www.doxygen.nl/manual/markdown.html)
