from collections import UserDict as DictClass
from pickle import dumps, loads, HIGHEST_PROTOCOL as PICKLE_PROTOCOL
import sqlite3
import os
import tempfile
from typing import Iterable, Tuple

def encode(obj):
    """Serialize an object using pickle to a binary format accepted by SQLite."""
    return sqlite3.Binary(dumps(obj, protocol=PICKLE_PROTOCOL))


def decode(obj):
    """Deserialize objects retrieved from SQLite."""
    return loads(bytes(obj))


class SqliteHandleDict(DictClass):
    VALID_FLAGS = ['c', 'r', 'w', 'n']
    
    def __init__(self, filename=None, tablename='unnamed', flag='c',
                autocommit = False, journal_mode = 'DELETE',
                encode=encode, decode=decode, timeout=5, outer_stack=True):
        """
        The `flag` parameter. Exactly one of:
          'c': default mode, open for read/write, creating the db/table if necessary.
          'w': open for r/w, but drop `tablename` contents first (start with empty table)
          'r': open as read-only
          'n': create a new database (erasing any existing tables, not just `tablename`!).

        The `encode` and `decode` parameters are used to customize how the values
        are serialized and deserialized.
        The `encode` parameter must be a function that takes a single Python
        object and returns a serialized representation.
        The `decode` function must be a function that takes the serialized
        representation produced by `encode` and returns a deserialized Python
        object.
        The default is to use pickle.

        The `timeout` defines the maximum time (in seconds) to wait for initial Thread startup.

        """
        self.in_temp = filename is None
        if self.in_temp:
            fd, filename = tempfile.mkstemp(prefix='sqldict')
            os.close(fd)

        if flag not in SqliteHandleDict.VALID_FLAGS:
            raise RuntimeError("Unrecognized flag: %s" % flag)
        self.flag = flag

        if flag == 'n':
            if os.path.exists(filename):
                os.remove(filename)

        dirname = os.path.dirname(filename)
        if dirname:
            if not os.path.exists(dirname):
                raise RuntimeError('Error! The directory does not exist, %s' % dirname)

        self.filename = filename

        # Use standard SQL escaping of double quote characters in identifiers, by doubling them.
        # See https://github.com/RaRe-Technologies/sqlitedict/pull/113
        self.tablename = tablename.replace('"', '""')

        self.autocommit = autocommit
        self.journal_mode = journal_mode
        self.encode = encode
        self.decode = decode
        self.timeout = timeout
        self._outer_stack = outer_stack

        connectionString = filename
        isolation_lvl = ''
        if self.flag == 'r':
            connectionString = 'file:' + connectionString + '?mode=ro'
            isolation_lvl = None
        

        self.conn: sqlite3.Connection = sqlite3.connect(connectionString, timeout,\
            isolation_level=isolation_lvl, check_same_thread=False)
        if self.flag == 'r':
            if self.tablename not in SqliteHandleDict.get_tablenames(self.filename):
                msg = 'Refusing to create a new table "%s" in read-only DB mode' % tablename
                raise RuntimeError(msg)
        else:
            MAKE_TABLE = 'CREATE TABLE IF NOT EXISTS "%s" (key TEXT PRIMARY KEY, value BLOB)' % self.tablename
            self.conn.execute(MAKE_TABLE)
            self.conn.commit()
        if flag == 'w':
            self.clear()
    
    @staticmethod
    def get_tablenames(filename):
        """get the names of the tables in an sqlite db as a list"""
        if not os.path.isfile(filename):
            raise IOError('file %s does not exist' % (filename))
        GET_TABLENAMES = 'SELECT name FROM sqlite_master WHERE type="table"'
        with sqlite3.connect(filename) as conn:
            cursor = conn.execute(GET_TABLENAMES)
            res = cursor.fetchall()

        return [name[0] for name in res]
    
    def commit(self, blocking=True):
        self.conn.commit()
    
    def __len__(self):
        COUNT_ALL = 'SELECT COUNT(*) FROM, %s' % self.filename
        count = self.conn.execute(COUNT_ALL).fetchone()
        return count[0] if count is not None else 0
            
    
    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        if self.autocommit:
            self.commit()
        
        self.close()
        
    def close(self):
        if self.autocommit:
            self.commit()
        self.conn.close()
    
    def clear(self):
        if self.flag == 'r':
            raise RuntimeError('Refusing to clear read-only SqliteDict')

        # avoid VACUUM, as it gives "OperationalError: database schema has changed"
        CLEAR_ALL = 'DELETE FROM "%s";' % self.tablename
        self.conn.commit()
        self.conn.execute(CLEAR_ALL)
        self.conn.commit()
    
    def terminate(self):
        if self.flag == 'r':
            raise RuntimeError('Refusing to terminate database from read-only connection')

        self.close()

        if self.filename == ':memory:':
            return

        try:
            if os.path.isfile(self.filename):
                os.remove(self.filename)
        except (OSError, IOError):
            print("failed to delete %s" % (self.filename))
        


    def set(self, key, item):
        if self.flag == 'r':
            raise RuntimeError('Refusing to terminate database from read-only connection')
        
        INSERT_ENTRY = 'INSERT INTO "%s" (key, value) VALUES (?,?)' % self.tablename
        self.conn.execute(INSERT_ENTRY, (key, self.encode(item)))
        if self.autocommit:
            self.commit()
    
    def update(self, items: Iterable[Tuple[any, any]]):
        INSERT_MANY_ENTRY = 'INSERT INTO "%s" (key, value) VALUES (?,?)' % self.tablename
        lst = ( (key, self.encode(value)) for key, value in items)
        self.conn.executemany(INSERT_MANY_ENTRY, lst)
        
        if self.autocommit:
            self.commit()

    DEFAULT_ERRORVALUE = '__ERRORVALUE__'
    def __getitem__(self, key):
        value = self.get(key, default_value=SqliteHandleDict.DEFAULT_ERRORVALUE)
        if value is SqliteHandleDict.DEFAULT_ERRORVALUE:
            raise KeyError(key)
        return value
    
    def get(self, key, default_value = None):
        GET_ITEM = 'SELECT value FROM "%s" WHERE key = ?' % self.tablename
        item = self.conn.execute(GET_ITEM, (key,)).fetchone()
        
        if item is None:
            return default_value
        
        return self.decode(item[0])
        

    def __setitem__(self, key, item):
        self.set(key, item)

    def insertOrReplace(self, key, item):
        INSERT_ENTRY = 'REPLACE INTO "%s" (key, value) VALUES (?,?)' % self.tablename
        self.conn.execute(INSERT_ENTRY, (key, self.encode(item)))
        if self.autocommit:
            self.commit()

    def __delitem__(self, key):
        if self.flag == 'r':
            raise RuntimeError('Refusing to delete from read-only SqliteDict')

        if key not in self:
            raise KeyError(key)
        DEL_ITEM = 'DELETE FROM "%s" WHERE key = ?' % self.tablename
        self.conn.execute(DEL_ITEM, (key,))
        if self.autocommit:
            self.commit()

    def __iter__(self):
        return self.iterkeys()
    
    def iterkeys(self):
        GET_KEYS = 'SELECT key FROM "%s" ORDER BY rowid' % self.tablename
        for key in self.conn.execute(GET_KEYS).fetchall():
            yield key[0]

    def itervalues(self):
        GET_VALUES = 'SELECT value FROM "%s" ORDER BY rowid' % self.tablename
        for value in self.conn.execute(GET_VALUES).fetchall():
            yield self.decode(value[0])

    def iteritems(self):
        GET_ITEMS = 'SELECT key, value FROM "%s" ORDER BY rowid' % self.tablename
        for key, value in self.conn.execute(GET_ITEMS).fetchall():
            yield key, self.decode(value)

    def keys(self):
        return self.iterkeys()

    def values(self):
        return self.itervalues()

    def items(self):
        return self.iteritems()

    def __contains__(self, key):
        HAS_ITEM = 'SELECT 1 FROM "%s" WHERE key = ?' % self.tablename
        return self.conn.execute(HAS_ITEM, (key,)).fetchone() is not None

    # Now, add the methods in dicts but not in MutableMapping
    def __repr__(self):
        return repr(self)

    def __or__(self, other):
        raise sqlite3.NotSupportedError("this function is not supported")

    def __ror__(self, other):
        raise sqlite3.NotSupportedError("this function is not supported")

    def __ior__(self, other):
        raise sqlite3.NotSupportedError("this function is not supported")

    def __copy__(self):
        raise sqlite3.NotSupportedError("this function is not supported")

    def copy(self):
        raise sqlite3.NotSupportedError("this function is not supported")
    