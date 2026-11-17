SCHEMA_VERSION = 1

TASKS_COLUMNS = """(
        task_id INTEGER PRIMARY KEY AUTOINCREMENT,
        product_id TEXT NOT NULL,
        revision INTEGER NOT NULL,
        digest TEXT NOT NULL,
        calculation TEXT NOT NULL CHECK(calculation IN
            ('primitive-v1','connectivity-v1')),
        FOREIGN KEY(product_id,revision)
            REFERENCES revisions(product_id,revision))"""

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
    "CREATE TABLE IF NOT EXISTS tasks " + TASKS_COLUMNS,
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


def migrate_tasks_table(connection):
    """Widen the V1 tasks CHECK in place, preserving every row.

    V1 databases keep user_version=1; only the tasks table definition is
    rebuilt so it also accepts 'connectivity-v1'. Called inside the
    initialization transaction with foreign keys temporarily off; the
    foreign_key_check afterwards guards the copy.
    """
    row = connection.execute(
        "SELECT sql FROM sqlite_master WHERE type='table' AND name='tasks'"
    ).fetchone()
    if row is None or "connectivity-v1" in row[0]:
        return
    connection.execute("CREATE TABLE tasks_migration " + TASKS_COLUMNS)
    connection.execute(
        "INSERT INTO tasks_migration"
        "(task_id,product_id,revision,digest,calculation) "
        "SELECT task_id,product_id,revision,digest,calculation FROM tasks")
    connection.execute("DROP TABLE tasks")
    connection.execute("ALTER TABLE tasks_migration RENAME TO tasks")


def initialize_schema(connection):
    check_schema_version(connection)
    connection.execute("PRAGMA foreign_keys=OFF")
    try:
        connection.execute("BEGIN IMMEDIATE")
        try:
            for statement in TABLES:
                connection.execute(statement)
            migrate_tasks_table(connection)
            connection.execute("PRAGMA user_version=1")
            if connection.execute("PRAGMA foreign_key_check").fetchall():
                raise ValueError("schema initialization broke foreign keys")
            connection.execute("COMMIT")
        except BaseException:
            if connection.in_transaction:
                connection.execute("ROLLBACK")
            raise
    finally:
        connection.execute("PRAGMA foreign_keys=ON")
