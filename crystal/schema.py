SCHEMA_VERSION = 1

TABLES = (
    """CREATE TABLE IF NOT EXISTS products (
        product_id TEXT PRIMARY KEY,
        current_revision INTEGER,
        FOREIGN KEY(product_id,current_revision)
            REFERENCES revisions(product_id,revision)
            DEFERRABLE INITIALLY DEFERRED)""",
    """CREATE TABLE IF NOT EXISTS revisions (
        product_id TEXT NOT NULL REFERENCES products(product_id),
        revision INTEGER NOT NULL CHECK(revision > 0),
        digest TEXT NOT NULL,
        canonical_json TEXT NOT NULL,
        raw_json TEXT NOT NULL,
        PRIMARY KEY(product_id,revision))""",
    """CREATE TABLE IF NOT EXISTS tasks (
        task_id INTEGER PRIMARY KEY AUTOINCREMENT,
        product_id TEXT NOT NULL,
        revision INTEGER NOT NULL,
        digest TEXT NOT NULL,
        calculation TEXT NOT NULL CHECK(calculation = 'primitive-v1'),
        FOREIGN KEY(product_id,revision)
            REFERENCES revisions(product_id,revision))""",
    """CREATE TABLE IF NOT EXISTS results (
        result_id INTEGER PRIMARY KEY AUTOINCREMENT,
        task_id INTEGER NOT NULL UNIQUE REFERENCES tasks(task_id))""",
    """CREATE TABLE IF NOT EXISTS result_bodies (
        result_id INTEGER PRIMARY KEY REFERENCES results(result_id),
        payload TEXT NOT NULL)""",
)


def check_schema_version(connection):
    version = connection.execute("PRAGMA user_version").fetchone()[0]
    if version not in (0, SCHEMA_VERSION):
        raise ValueError("incompatible database schema")
    if version == 0:
        tables = {row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%'")}
        if tables:
            raise ValueError("unversioned nonempty database")
    return version


def initialize_schema(connection):
    check_schema_version(connection)
    connection.execute("BEGIN IMMEDIATE")
    try:
        for statement in TABLES:
            connection.execute(statement)
        connection.execute("PRAGMA user_version=1")
        connection.execute("COMMIT")
    except BaseException:
        if connection.in_transaction:
            connection.execute("ROLLBACK")
        raise
