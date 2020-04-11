import os

from sqlalchemy import create_engine

DATABASE_HOST = os.environ.get('DATABASE_HOST', 'localhost')
DATABASE_PORT = os.environ.get('DATABASE_PORT', '5432')

RAZI_DB = os.environ.get('RAZI_DB', 'razi_test')
RAZI_USER = os.environ.get('RAZI_USER', 'user')
RAZI_PASSWORD = os.environ.get('RAZI_PASSWORD', 'password')

engine = create_engine(f'postgresql://{RAZI_USER}:{RAZI_PASSWORD}@{DATABASE_HOST}:{DATABASE_PORT}/{RAZI_DB}')
