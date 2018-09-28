"""
    Module to facilitate the database access and some common methods used throughout code
"""


from sqlalchemy.sql import select
from sqlalchemy import Table
from sqlalchemy import (MetaData, create_engine)
from sqlalchemy.engine import reflection
from sqlalchemy.schema import CreateSchema

from src import sqlalchemy_extension as se
import utils


class DbConnection(metaclass=utils.Singleton):
    """
    This class is responsible to facilitate the connection and interaction with the database.
    """
    def __init__(self):
        self.eng = None
        self.meta = None
        self.schema_output = None

        self._is_connection_initialized = False

    def db_init(self, url_connection, schema_output=None):
        """
            Args:
                url_connection: An object returned by the class sqlalchemy.engine.url.URL
                schema_output: schema to save the tables.
        """
        self.eng = create_engine(url_connection)
        self.meta = MetaData(bind=self.eng)
        self.schema_output = schema_output

        if self.schema_output:
            self._create_schema_if_not_exist()

        self._is_connection_initialized = True

    def is_connection_initialized(self):
        return self._is_connection_initialized

    def select_columns(self, table_name, schema=None, columns=None):
        """
            It returns a tuple of values defined by a list of columns. If the schema is not defined, the schema_output
        is used.
        """
        if schema is None:
            schema = self.schema_output
        with self.eng.connect() as con:
            table = Table(table_name, MetaData(bind=con),
                          autoload=True, schema=schema)
            cols = self.create_columns_sql_format(table, columns)
            stm = select(cols).select_from(table)
            result = con.execute(stm)
            return result.fetchall()

    def prepare_data_to_plot_catalog(self, table_dependencies):
        """
            The catalogs need the ra, dec info to generate plots. When the table doesn't have this info, it must look
        for this information on the dataset.
        """
        with self.eng.connect() as con:
            t_table_name = Table(table_dependencies['cur_table']['table'], MetaData(bind=con),
                                 autoload=True, schema=table_dependencies['cur_table']['schema'])
            t_dataset = Table(table_dependencies['coadd_objects']['table'], MetaData(bind=con),
                              autoload=True, schema=table_dependencies['coadd_objects']['schema'])

            stm_join = t_dataset
            stm_join = stm_join.join(
                t_table_name, t_dataset.c.coadd_objects_id == t_table_name.c.coadd_objects_id)

            stm = select([t_dataset.c.ra, t_dataset.c.dec]).select_from(stm_join)
            result = con.execute(stm)
            return result.fetchall()

    def _create_schema_if_not_exist(self):
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
        # runs a transaction
        with self.eng.begin() as con:
            con.execute(se.CreateTableAs(self.schema_output, table_name, stm))

    def delete_table(self, table_name):
        with self.eng.begin() as con:
            con.execute(se.DropTable(self.schema_output, table_name))
