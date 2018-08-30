from sqlalchemy import Integer, func
from sqlalchemy.sql import select
from sqlalchemy.sql.expression import bindparam, ClauseElement, Executable
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.sql.elements import ColumnElement, Visitable, Grouping
from sqlalchemy.sql.operators import custom_op, is_precedent


"""
    Some Sql operations are written differently based on the backend used. So,
the classes below give support for this cases.
"""


class BitwiseAnd(ColumnElement):
    """
    Bitwise and operation.
    """
    type = Integer()
    operator = custom_op("&", precedence=6)

    def __init__(self, left, right):
        if not isinstance(left, Visitable):
            left = bindparam("bitand", left, unique=True)
        if not isinstance(right, Visitable):
            right = bindparam("bitand", right, unique=True)
        self.left = left
        self.right = right

    def self_group(self, against=None):
        if is_precedent(self.operator, against):
            return Grouping(self)
        else:
            return self


@compiles(BitwiseAnd)
def _compile_bitwise_and(element, compiler, **kwargs):
    left = element.left.self_group(against=element.operator)
    right = element.right.self_group(against=element.operator)
    return compiler.process(element.operator(left, right))


@compiles(BitwiseAnd, "oracle")
def _compile_bitwise_and_oracle(element, compiler, **kwargs):
    return compiler.process(func.BITAND(element.left, element.right))


class CreateTableAs(Executable, ClauseElement):
    """
    Creates a new table in the database using a query result.
    """
    def __init__(self, schema, name, query):
        self.schema = schema
        self.name = name
        self.query = query


@compiles(CreateTableAs)
def _create_table_as(element, compiler, **kw):
    _schema = "%s." % element.schema if element.schema is not None else ''
    return "CREATE TABLE %s%s AS (%s)" % (
        _schema,
        element.name,
        compiler.process(element.query)
    )


class DropTable(Executable, ClauseElement):
    """
    Drop a table in the database.
    """
    def __init__(self, schema, name):
        self.schema = schema
        self.name = name


@compiles(DropTable)
def _drop_table(element, compiler, **kw):
    _schema = "%s." % element.schema if element.schema is not None else ''
    return "DROP TABLE %s%s" % (_schema, element.name)


# @compiles(DropTable, "oracle")
# def _drop_table(element, compiler, **kw):
#     return "DROP TABLE %s PURGE" % (element.name)


def create_columns_sql_format(table_obj, columns):
    t_columns = table_obj
    if columns is not None:
        t_columns = list()
        for col in columns:
            t_columns.append(getattr(table_obj.c, col))
    return t_columns


def select_columns(engine, table, columns=None):
    with engine.connect() as con:
        cols = create_columns_sql_format(table, columns)
        stm = select(cols).select_from(table)
        result = con.execute(stm)
        return result.fetchall()
