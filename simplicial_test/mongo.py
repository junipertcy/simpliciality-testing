from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError, NetworkTimeout
import configparser
import hashlib


class DB(object):
    def __init__(self, db, col, degs, sizes):
        config = configparser.ConfigParser()
        config.read('/Users/tzuchi/.config/mongoconfig/master.ini')
        m_config = config['mongo.netscied.tw']

        client = MongoClient(
            'mongodb://admin:' + m_config['pwd'] + '@' + m_config['ip'] + ':' + m_config['port'] + '/?authSource=admin')
        self.col = client[db][col]
        total = int(sum(degs))
        m = len(sizes)
        n = len(degs)
        self.total = total
        self.m = m
        self.n = n
        self.degs = degs
        self.sizes = sizes
        str2hash = str(degs) + str(sizes)
        _hash = hashlib.sha3_512(str2hash.encode())
        self.hash = _hash.hexdigest()[-8:]
        self.key = {
            "total": self.total,
            "m": self.m,
            "n": self.n,
            "degs": self.degs,
            "sizes": self.sizes,
            "sha3_512_8": self.hash
        }

    def init(self):
        try:
            self.col.find_one_and_update(self.key, {'$set': {'count': 0, 'simplicial': False}}, upsert=True)
        except NetworkTimeout:
            pass
        except DuplicateKeyError:
            pass

    def add_to_failed_attempts(self):
        self.col.find_one_and_update({
            "total": self.total,
            "sha3_512_8": self.hash
        }, {'$inc': {'count': 1}, '$set': {'simplicial': False}}, upsert=True)

    def mark_simplicial(self):
        self.col.find_one_and_update({
            "total": self.total,
            "sha3_512_8": self.hash
        }, {'$set': {'simplicial': True}}, upsert=True)
