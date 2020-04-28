#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $script_dir
jupyter nbconvert --to markdown plot_chain.ipynb
# replace files relative location
sed 's#plot_chain_files#scripts\/notebooks\/plot_chain_files#g' -i plot_chain.md

source_dir="$script_dir/../../"
cd $source_dir
# Delete from Notebook example until the end of the file
sed '/Notebook example/,$d' -i README.md
# # Append to current README.md
echo -e "## Notebook example\n" >> README.md
cat ./scripts/notebooks/plot_chain.md >> README.md
echo "Do not append anything after Notebook header." >> README.md
# echo "Exit!"
