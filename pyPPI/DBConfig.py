"""Configuration file for access database
"""
import sys
from getpass import getpass

try:
    import MySQLdb
except:
    import pymysql as MySQLdb

__initlized = False

if sys.version_info[:2] <= (2, 7):
    input = raw_input

USER = ''
PASSWD = ''
DB_NAME = ''


def init_connection(user=None, passwd=None, db=None):
    """
    Initializes connection details to database
    :param user: database user
    :param passwd: database password
    :param db: database name
    """
    global USER, PASSWD, DB_NAME, __initlized
    if user is None:
        user = input('Enter db user: ')
    if passwd is None:
        passwd = getpass('password: ')
    if db is None:
        db = input('Enter db name: ')
    USER = user
    PASSWD = passwd
    DB_NAME = db
    __initlized = True


def get_connection():
    """
    Get connection object to local database
    """
    global USER, PASSWD, DB_NAME, __initlized
    if not __initlized:
        init_connection()
    return MySQLdb.connect(host='localhost', user=USER, passwd=PASSWD, db=DB_NAME, local_infile=1)
