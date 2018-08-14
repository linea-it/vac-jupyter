"""
    Access to the database, metadata and some methods.
"""


from sqlalchemy.sql import select
from sqlalchemy import Table
from sqlalchemy import (MetaData, create_engine)
from sqlalchemy.engine import reflection
from sqlalchemy.schema import CreateSchema

from src import sqlalchemy_extension as se


class DataAccessLayer():
    def __init__(self, url, schema_output):
        self.eng = create_engine(url)
        self.meta = MetaData(bind=self.eng)
        self.schema_output = schema_output

        self.create_schema_if_not_exist()

    def select_columns(self, table_name, schema=None, columns=None):
        if schema is None:
            schema = self.schema_output
        with self.eng.connect() as con:
            table = Table(table_name, self.meta,
                          autoload=True, schema=schema)
            cols = self.create_columns_sql_format(table, columns)
            stm = select(cols).select_from(table)
            result = con.execute(stm)
            return result.fetchall()

    def prepare_data_to_plot_catalog(self, table_dependencies):
        with self.eng.connect() as con:
            t_table_name = Table(table_dependencies['cur_table']['table'], self.meta,
                                 autoload=True, schema=table_dependencies['cur_table']['schema'])
            t_dataset = Table(table_dependencies['coadd_objects']['table'], self.meta,
                              autoload=True, schema=table_dependencies['coadd_objects']['schema'])

            stm_join = t_dataset
            stm_join = stm_join.join(
                t_table_name, t_dataset.c.coadd_objects_id == t_table_name.c.coadd_objects_id)

            stm = select([t_dataset.c.ra, t_dataset.c.dec]).select_from(stm_join)
            result = con.execute(stm)
            return result.fetchall()

    def create_schema_if_not_exist(self):
        insp = reflection.Inspector.from_engine(self.eng)
        if self.schema_output not in insp.get_schema_names():
            self.eng.execute(CreateSchema(self.schema_output))

    @staticmethod
    def create_columns_sql_format(table_obj, columns):
        t_columns = table_obj
        if columns is not None:
            t_columns = list()
            for col in columns:
                t_columns.append(getattr(table_obj.c, col))
        return t_columns

    def create_table(self, table_name, stm):
        with self.eng.connect() as con:
            con.execute("commit")
            con.execute(se.CreateTableAs(self.schema_output, table_name, stm))

    def delete_table(self, table_name):
        with self.eng.connect() as con:
            con.execute("commit")
            con.execute(se.DropTable(self.schema_output, table_name))
