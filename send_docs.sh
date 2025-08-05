python clear_comments.py
rm docs/code/*
cp clean/* docs/code/

# Plotting
cp data_convection/anim/*.gif docs/gif/
cp plot_*.py docs/scripts/

cd docs/code
zip MilleFEuiIle *.py
cd ..

make clean
make html
cp -r _build/html/* /nfs1/WWW/users/kihoulou/millefeuille
cp code/MilleFEuiIle.zip /nfs1/WWW/users/kihoulou/millefeuille