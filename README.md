1. Create virtual env
virtualenv --no-site-packages --always-copy --python python3 env

pip install sqlalchemy
pip install sqlparse
pip install psycopg2
pip install matplotlib
pip install healpy

git submodule init
git submodule sync
git submodule update
