#/usr/bin/env python3
#--*--coding:utf-8--*--
import pandas as pd
import mysql.connector
from sqlalchemy import create_engine
from sqlalchemy.types import VARCHAR
import numpy as np
df = pd.read_csv('./sed008_specimen.csv',
                 skiprows=0,
                 header=1,
                 delimiter=',',
                 index_col=False)
print(df)
df.replace(999,np.NaN)
#be care of the format of the dataframe is correct!!!!!!
engine = create_engine("mysql+mysqlconnector://root:root@127.0.0.1/geomagia", echo = False)
df.to_sql(name='geomagiaSed',con=engine, if_exists='replace',
          index=False,
          dtype={'None':VARCHAR(5)})
