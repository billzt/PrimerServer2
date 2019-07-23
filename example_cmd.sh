# design primers and check specificity
primertool full tests/query_design_multiple tests/example.fa -o full.json -t full.tsv

# design primers only
primertool design tests/query_design_multiple tests/example.fa -o design.json -t design.tsv

# check specificity only
primertool check tests/query_check_multiple tests/example.fa -o check.json -t check.tsv

