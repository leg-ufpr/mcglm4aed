# Removes the _site/.
rm -r _site/
# Rscript -e "library(rmarkdown); clean_site()"

# The root file.
echo $PWD
THEROOT=$PWD
cd $THEROOT

# List of Rmd files.
RMDFILES=$(find . -name \*.Rmd -print)
echo $RMDFILES

# Runs rmarkdown::render() in each Rmd file down in the tree.
for RMD in $RMDFILES; do
    DIRNAME=`dirname "$RMD"`
    RMDNAME=`basename "$RMD"`
    cd $DIRNAME
    Rscript -e "library(rmarkdown); render(\"$RMDNAME\")"
    cd $THEROOT
done

# Renders the site.
Rscript -e "library(rmarkdown); render_site()"

# Change to the gh-pages branch.
# git branch --delete gh-pages
git checkout -b gh-pages
# git checkout gh-pages

# Go to the _site/.
# cd _site/
# tree

# Force add html, js, css, png, jpg and site_libs/.
git status
git branch
git add -f site_libs/
git add -f slides/*.{html,css}
git add -f tutorials/*.{html,css}

git rm tutorials/*.Rmd
git rm slides/*.Rmd

# Undo the `git add`.
# git reset HEAD

git status
git commit -m "Updates gh-pages"

# git clean --dry-run -fdx
# git status

# Kill the gh-pages in the server.
git push origin --delete gh-pages

# Push the new gh-pages.
git push -u origin gh-pages

git checkout master
