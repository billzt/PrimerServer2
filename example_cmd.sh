# design primers and check specificity (the default mode)
primertool tests/query_design_multiple tests/example.fa -p 2 -o example_design_check.json -t example_design_check.tsv

# check specificity only
primertool tests/query_check_multiple tests/example.fa --only-specificity -p 2 -o example_check.json -t example_check.tsv
