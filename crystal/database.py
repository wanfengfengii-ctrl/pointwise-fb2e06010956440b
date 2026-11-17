import sqlite3
from contextlib import contextmanager
from .schema import initialize_schema


def open_database(path):
    connection = sqlite3.connect(path, isolation_level=None)
    try:
        connection.execute("PRAGMA foreign_keys=ON")
        initialize_schema(connection)
    except BaseException:
        connection.close()
        raise
    return connection


@contextmanager
def transaction(connection):
    if connection.in_transaction:
        raise RuntimeError("nested publication transaction is unsupported")
    connection.execute("BEGIN IMMEDIATE")
    try:
        yield connection
        connection.execute("COMMIT")
    except BaseException:
        if connection.in_transaction:
            connection.execute("ROLLBACK")
        raise


def close_database(connection):
    connection.close()
