# Shared defaults for command-line and editor builds run from this directory.
$pdf_mode = 1;
$out_dir = 'build';
# TeX cannot create the subdirectory used by \include's auxiliary files.
use File::Path qw(make_path);
make_path('build/chapters');
# These switches also survive latexmk's -usepretex command replacement.
$pdflatex_default_switches = '-interaction=nonstopmode -halt-on-error -file-line-error -synctex=1';
